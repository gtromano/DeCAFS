#' Estimate parameter in the Random Walk Autoregressive model
#' 
#' 
#' This function perform robust estimation of parameters in the Random Walk plus Autoregressive model using a method of moments estimator. To model the time-dependency DeCAFS relies on three parameters. These are \code{sdEta}, the standard deviation of the drift (random fluctuations) in the signal, modeled as a Random Walk process, \code{sdNu}, the standard deviation of the AR(1) noise process, and \code{phi}, the autocorrelation parameter of the noise process. 
#' The final estimation of the change locations is affected by the l0 penalty beta and the estimation of the process by those three initial parameters. Therefore, the choice of penalties for DeCAFS is important: where possible investigate resulting segmentations. Should the algorithm return a misspecified estimation of the signal, it might be good to constrain the estimation of the parameters to an edge case. This can be done through the argument \code{model}. Alternatively, one could employ a range of penalties or tune these on training data. To manually specify different penalties, see \code{\link[=DeCAFS]{DeCAFS()}} documentation.
#' If unsure of which model is the most suited for a given sequence, see  \code{\link[=guidedModelSelection]{guidedModelSelection()}} for guided model selection.
#'
#' @param y A vector of observations
#' @param model Constrain estimation to an edge case of the RWAR model. Defaults to \code{"RWAR"}. To fit an AR model only with a piece-wise constant signal, specify \code{"AR"}. To fit a a random walk plus noise, specify \code{"RW"}.
#' @param K The number of K-lags differences of the data to run the robust estimation over. Default set at 15.
#' @param phiLower Smallest value of the autocorrelation parameter. Default set at 0.
#' @param phiUpper Highest value of the autocorrelation parameter. Default set at 0.99.
#' @param sdEtaUpper Highest value of the RW standard deviation. Default set at Inf
#' @param sdNuUpper Highest value of the AR(1) noise standard deviation. Default set at Inf
#' @param warningMessage A message to warn the user when the automatic parameter estimation is employed.
#'
#' @return 
#'  Returns a list of estimates that can be employed as an argument for parameter \code{modelParam} to run \code{\link[=DeCAFS]{DeCAFS()}}. Those are:
#' \describe{
#' \item{\code{sdEta}}{the SD of the drift (random fluctuations) in the signal,}
#' \item{\code{sdNu}}{the SD of the AR(1) noise process,}
#' \item{\code{phi}}{the autocorrelation parameter of the noise process.}
#' }
#' 
#' @export
#'
#' @examples
#' set.seed(42)
#' y <- dataRWAR(n = 1e3, phi = .5, sdEta = 1, sdNu = 3,  jumpSize = 15, type = "updown", nbSeg = 5)$y
#' estimateParameters(y)
estimateParameters <- function (y, model = c("RWAR", "AR", "RW"), K = 15, phiLower = 0, phiUpper = .999, sdEtaUpper = Inf, sdNuUpper = Inf, warningMessage = FALSE) 
{
  
  if (warningMessage) warning("\nAutomatic parameters estimation employed.\nThe choice of penalties for DeCAFS is important - where possible investigate the returned segmentation.\nFor more details about initial parameters estimation, please see help('estimateParameters').\nIf unsure of which model is the most suited for a given sequence, run guidedModelSelection().\n\nTo silence this warning when using DeCAFS, run DeCAFS with argument: warningMessage = FALSE.", call. = F)
  
  model <- match.arg(model)
  
  if (model == "RWAR") {
    n <- length(y)
    if (!is.numeric(y)) 
      stop("Please provvide a vector of observations y")
    if (K > (n + 1)) {
      K <- n - 1
      warning(paste0("Lag parameter K is too big. Setting lag to n-1, i.e.: ", K))
    }
    
    # initial estimates for our sigmaEta, sigmaNu
    phi <- 0.5
    Wk2 <- sapply(1:K, function(k) {
      zk <- y[(1 + k):n] - y[1:(n - k)]
      mad(zk)^2
    })
    nuX <- sapply(1:K, function(k) {
      2 * (1 - phi^k)/(1 - phi^2)
    })
    etaX <- 1:K
    model <- lm(Wk2 ~ -1 + etaX + nuX)
    sdEta <- abs(coefficients(model)[1])
    sdNu <- abs(coefficients(model)[2])
    
    start = c(max(0, phiLower), sdNu, sdEta)
    out = optim(
      par = start,
      .MoMCost,
      lower = c(phiLower, 0.001, 0),
      upper = c(phiUpper, sdNuUpper, sdEtaUpper),
      V = Wk2,
      method = "L-BFGS-B"
    )
    
    return(list(sdEta = as.numeric(sqrt(out$par[3])), sdNu = as.numeric(sqrt(out$par[2])), phi =  out$par[1]))
  } else if (model == "AR") {
    est <- estimateParameters(y, sdEtaUpper = 1e-12)
    est$sdEta <- 0
    return(est)
  } else {
    est <- estimateParameters(y, phiUpper =  1e-12)
    est$phi <- 0
    return(est)
  }
  
}

.MoMCost <- function(par, V) {
  means <- sapply(1:length(V), function(k)
    k * par[3] + 2 * par[2] * (1 - par[1] ^ k)/(1 - par[1]^2))
  return(sum((means - V) ^ 2))
}