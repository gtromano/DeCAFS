#' Estimate parameter in the Random Walk Autoregressive model
#' 
#' This function perform robust estimation of parameters in the Random Walk plus Autoregressive model using
#' a method of moments estimator. 
#' Returns a list of estimates that can be employed as an argument for parameter \code{modelParam} to run \code{\link[=DeCAFS]{DeCAFS()}}.
#' 
#' @param y A vector of observations
#' @param K The number of K-lags differences of the data to run the robust estimation over. Default set at 15.
#' @param phiLower Smallest value of the autocorrelation parameter. Default set at 0.
#' @param phiUpper Highest value of the autocorrelation parameter. Default set at 0.99.
#' @param sdEtaUpper Highest value of the RW standard deviation. Default set at Inf
#' @param sdNuUpper Highest value of the AR(1) noise standard deviation. Default set at Inf
#' @param model Constrain estimation to an edge case of the RWAR model. Defaults to \code{"RWAR"}. To fit an AR model only with a piece-wise constant signal, specify \code{"AR"}. To fit a a random walk plus noise, specify \code{"RW"}.
#' 
#' @return 
#' A list containing:
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

estimateParameters <- function (y, K = 15, phiLower = 0, phiUpper = .999, sdEtaUpper = Inf, sdNuUpper = Inf, model = c("RWAR", "AR", "RW")) 
{
  
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