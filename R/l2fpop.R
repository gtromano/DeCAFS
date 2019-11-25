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