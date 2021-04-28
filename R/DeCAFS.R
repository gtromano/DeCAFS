#' Main DeCAFS algorithm for detecting abrupt changes
#' 
#' This function implements the DeCAFS algorithm to detect abrupt changes in mean of a univariate data stream in the presence of local fluctuations and auto-correlated noise. 
#' It detects the changes under a penalised likelihood model where the data, \eqn{y_1, ..., y_n}, is \deqn{y_t = \mu_t + \epsilon_t}
#' with \eqn{\epsilon_t} an AR(1) process, and for \eqn{t = 2, ..., N} 
#' \deqn{\mu_t = \mu_{t-1} + \eta_t + \delta_t}
#' where at time \eqn{t} if we do not have a change then \eqn{\delta_t = 0} and \eqn{\eta_t \sim N(0, \sigma_\eta^2)}; whereas if we have a change then \eqn{\delta_t \neq 0} and \eqn{\eta_t=0}.
#' DeCAFS estimates the change by minimising a cost equal to twice the negative log-likelihood of this model, with a penalty \eqn{\beta} for adding a change.
#' Note that the default DeCAFS behavior will assume the RWAR model, but fit on edge cases is still possible. For instance, should the user wish for DeCAFS to fit an AR model only with a piecewise constant signal, or similarly a model that just assumes random fluctuations in the signal, this can be specified within the initial parameter estimation, by setting the argument: \code{modelParam = estimateParameters(y, model = "AR")}. Similarly, to allow for negative autocorrelation estimation, set \code{modelParam = estimateParameters(Y$y, phiLower = -1)}.
#' 
#' 
#' @param data A vector of observations y
#' @param beta The l0 penalty. The default one is \code{2 * log(N)} where \code{N} is the length of the data.
#' @param modelParam A list of 3 initial model parameters: \code{sdEta}, the SD of the drift (random fluctuations) in the signal, \code{sdNu}, the SD of the AR(1) noise process, and \code{phi}, the autocorrelation parameter of the noise process (so the stationary variance of the AR(1) noise process is \code{sdnu^2} / (1 - \code{phi^2}). Defaulted to \code{estimateParameters(data, K = 15)}, to perform automatically estimation of the three. See \code{\link[=estimateParameters]{estimateParameters()}} for more details.
#' @param penalties Can be used as an alternative to the model parameters, a list of 3 initial penalties: \code{lambda}, the l2-penalty penalising over the lag-1 of the signal, \code{gamma}, penalising over the lag-1 of the AR(1) noise process, \code{phi}, the autocorrelation parameter. These are related to the \code{modelParam} list by \code{list(lambda = 1 / sdEta ^ 2, gamma = 1 / sdNu ^ 2, phi = phi)}. Only one argument between \code{penalties} and \code{modelParam} should be specified. Defaulted to NULL. 
#' @param type The type of change one wants to look for. At the moment only 'std' is implemented.
#'
#' @return Returns an s3 object of class DeCAFSout where:
#' \describe{
#' \item{\code{$changepoints}}{is the vector of change-point locations,}
#' \item{\code{$signal}}{is the estimated signal without the auto-correlated noise,}
#' \item{\code{$costFunction}}{ is the optimal cost in form of piecewise quadratics at the end of the sequence,} 
#' \item{\code{$estimatedParameters}}{is a list of parameters estimates (if estimated, otherwise simply the initial \code{modelParam} input),} 
#' \item{\code{$data}}{is the sequence of observations.}
#' }
#' 
#' @export
#'
#' @references Romano, G., Rigaill, G., Runge, V., Fearnhead, P. (2021). Detecting Abrupt Changes in the Presence of Local Fluctuations and Autocorrelated Noise. Journal of the American Statistical Association. \doi{10.1080/01621459.2021.1909598}.
#'
#' @examples
#' library(ggplot2)
#' set.seed(42)
#' Y <- dataRWAR(n = 1e3, phi = .5, sdEta = 1, sdNu = 3, jumpSize = 15, type = "updown", nbSeg = 5)
#' y <- Y$y
#' res = DeCAFS(y)
#' ggplot(data.frame(t = 1:length(y), y), aes(x = t, y = y)) +
#'   geom_point() +
#'   geom_vline(xintercept = res$changepoints, color = "red") +
#'   geom_vline(xintercept = Y$changepoints, col = "blue",  lty = 3)



DeCAFS <- function(data, beta = 2 * log(length(data)), modelParam = estimateParameters(data), penalties = NULL, type = "std") {
  
  if(!is.numeric(data)) stop("Please provvide a vector of observations y")
  
  if (is.null(penalties)) {
    
    # here we use the model parameters (estimated by default trough estimateParameters)
    if(!is.list(modelParam)) stop("Please provvide a list of model parameters sdEta, sdNu, phi")
    pCheck <- c("sdEta", "sdNu", "phi") %in% names(modelParam)
    if (sum(pCheck) < 3) stop("Please provvide parameters: ", c("sdEta ", "sdNu ", "phi ")[!pCheck])
    
    lambda <- 1 / modelParam$sdEta ^ 2
    gamma <- 1 / modelParam$sdNu ^ 2
    phi <- modelParam$phi
    
  } else {
    
    if(!is.list(penalties)) stop("Please provvide a list of model penalties lambda, gamma, and autocorrelation parameter phi")
    pCheck <- c("lambda", "gamma", "phi") %in% names(penalties)
    if (sum(pCheck) < 3) stop("Please provvide penalties: ", c("lambda ", "gamma ", "phi ")[!pCheck])
    
    lambda <- penalties$lambda
    gamma <- penalties$gamma
    phi <- penalties$phi
    
  }
  
  
  # running the algorithm
  DeCAFSRes <- .DeCAFS(data, beta, lambda, gamma, phi, "std")
  
  # DeCAFSRes$changepoints <- DeCAFSRes$changepoints[-length(DeCAFSRes$changepoints)] - 1
  
  output <- list(changepoints = DeCAFSRes$changepoints,
                 signal = DeCAFSRes$signal,
                 costFunction = DeCAFSRes$costFunction,
                 modelParameters = modelParam,
                 data = data)
  class(output) <-  c("DeCAFSout", class(output))
  return(output)
}
