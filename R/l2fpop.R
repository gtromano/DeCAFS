#' Main l2fpop function
#' 
#' Detecting abrupt changes in mean in presence of autocorrelation or random fluctuations.
#' 
#' @param vectData A vector of observations y
#' @param beta The l0 penalty. The default one is 2 * log(N) where N is the lenght of the data.
#' @param lambda The l2 penalty for penalising over random fluctuations. When NULL it is estimated robustly from the data.
#' @param gamma The l2 penalty for penalising over autocorrelated noise. When NULL it is estimated robustly from the data.
#' @param phi The autocorrelation parameter phi. When NULL it is estimated robustly from the data. If one wish to exclude the autocorrelated component, then set phi = 0. In this case, then the gamma paramer is penalising over additional i.i.d. gaussian noise in the data.
#' @param type The type of change one wants to look for. At the moment only 'std' is implemented.
#'
#' @return
#' Returns a list where $changepoints is the vector of change-point locations, $costFunction is the optimal cost in form of piecewise quadratics at the end of the sequence, $estimatedParameters is a list of parameters estimates (if estimated), $data is the sequence of observations.
#' @export
#'
#' @examples
#' library(ggplot2)
#' set.seed(42)
#' Y = dataRWAR(n = 1e3, poisParam = .01, meanGap = 15, phi = .5, sdEta = 3, sdNu = 1)
#' y = Y$y
#' res = l2fpop(y)
#' ggplot(data.frame(t = 1:length(y), y), aes(x = t, y = y)) +
#'   geom_point() +
#'   geom_vline(xintercept = res$changepoints, color = 2) +
#'   geom_vline(xintercept = Y$changepoints, col = 4,  lty = 3)




l2fpop <- function(vectData, beta = 2 * log(length(vectData)), lambda = NULL, gamma = NULL, phi = NULL, type = "std") {
  
  if(!is.numeric(y)) stop("Please provvide a vector of observations y")
  
  estim <- NULL
  # estimating the parameter if they are not passed by the user
  if(length(c(lambda, gamma, phi)) != 3) {
    estim <- estimateParameters(y)
    if (is.null(lambda)) lambda <- 1 / estim$sdEta ^ 2
    if (is.null(gamma)) gamma <- 1 / estim$sdNu ^ 2
    if (is.null(phi)) phi <- estim$phi
  }
  
  # running the algorithm
  l2fpopRes <- .l2fpop(vectData, beta, lambda, gamma, phi, type)
  return(list(changepoints = l2fpopRes$changepoints,
              costFunction = l2fpopRes$costFunction,
              estimatedParameters = estim,
              data = y))
}