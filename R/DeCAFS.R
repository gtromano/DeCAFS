#' Main DeCAFS function
#' 
#' Detecting Abrupt Changes in the Presence ofLocal Fluctuations in the signal or Autocorrelated Noise.
#' 
#' @param data A vector of observations y
#' @param beta The l0 penalty. The default one is \code{2 * log(N)} where \code{N} is the length of the data.
#' @param modelParam A list of 3 initial model parameters: \code{sdEta}, the SD of the drift (random fluctuations) in the signal, \code{sdNu}, the SD of the AR(1) noise process, and \code{phi}, the autocorrelation parameter of the noise process. Defaulted to \code{estimateParameters(data, K = 15)}, to perform automatically estimation of the three. See \code{\link[=estimateParameters]{estimateParameters()}} for more details.
#' @param penalties Can be used as an alternative to the model parameters, a list of 3 initial penalties: \code{lambda} l2-penalty penalising over the lag-1 of the signal, \code{gamma}, penalising over the lag-1 of the AR(1) noise process, \code{phi}, the autocorrelation parameter. Defaulted to NULL. 
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
#' @references Romano, G., Rigaill, G., Runge, V., Fearnhead, P. Detecting Abrupt Changes in the Presence of Local Fluctuations and Autocorrelated Noise. arXiv preprint \url{https://arxiv.org/abs/2005.01379} (2020).
#'
#' @examples
#' library(ggplot2)
#' set.seed(42)
#' Y = dataRWAR(n = 1e3, poisParam = .01, meanGap = 15, phi = .5, sdEta = 1, sdNu = 3)
#' y = Y$y
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
  
  DeCAFSRes$changepoints <- DeCAFSRes$changepoints[-length(DeCAFSRes$changepoints)] - 1
  
  output <- list(changepoints = DeCAFSRes$changepoints,
                 signal = DeCAFSRes$signal,
                 costFunction = DeCAFSRes$costFunction,
                 modelParameters = modelParam,
                 data = data)
  class(output) <-  c("DeCAFSout", class(output))
  return(output)
}