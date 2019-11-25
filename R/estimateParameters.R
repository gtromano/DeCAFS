#' Estimate parameter in the Random Walk Autoregressive model
#' 
#' This function perform robust estimation of parameters in the Random Walk plus Autoregressive model using
#' a method of moments estimator.
#' 
#' @param y A vector of observations
#' @param K The number of lags to run the estimation over. Default set at 20. 
#'
#' @return 
#' A list containing $sigmaEta, the sd of the Random Walk Component, $sigmaNu, the sd of the AR noise, $phi, the autocorrelation parameter.
#' @export
#'
#' @examples
#' set.seed(42)
#' y <- dataRWAR(n = 1e4, poisParam = .01, phi = .7, sdEta = 4, sdNu = 3)$y
#' estimateParameters(y)

estimateParameters = function(y, K = 20) {
  n   <- length(y)
  
  if(!is.numeric(y)) stop("Please provvide a vector of observations y")
    
  if(K > (n+1)) {
    K <- n - 1
    warning(paste0("Lag parameter K is too big. Setting lag to n-1, i.e.: ", K))
  }
    
  phi <- .5
  
  # estimating the variances
  Wk2 <- sapply(1:K, function(k) {
    zk <- y[(1 + k):n] - y[1:(n-k)]
    zk <- zk[!zk %in% boxplot.stats(zk)$out] #  removing observations outside the 1.5 * IQR
    mean(zk^2)
  })
  
  nuX <- sapply(1:K, function(k) {
    2 * (1 - phi ^ k)/(1 - phi ^ 2)
  })
  
  etaX <- 1:K
  
  model <- lm(Wk2 ~ -1 + etaX + nuX)
  sdEta <- abs(coefficients(model)[1])
  sdNu <- abs(coefficients(model)[2])
  
  # re-estimating the phi
  p <- (Wk2[1] - sdEta)/(2 * sdNu)
  z <- y[(2):n] - y[1:(n-1)]
  phi <- c(-((p - 1)/p), # this is the solution to the equation of first grade
           mean(z[3:length(z)] * z[1:(length(z)-2)]) / mean(z[-1] * z[-length(z)]))
  phi <- mean(phi)
  
  if (is.nan(phi)) phi <- 0
  if (phi < 0 || phi > 1) {
    warning("Cannot estimate consistently autocorrelation parameter phi. Perharps no AR component is present in the data. Returning phi as 0.")
    phi <- 0
  }
  
  return(list(sdEta = as.numeric(sqrt(sdEta)), sdNu = as.numeric(sqrt(sdNu)), phi = as.numeric(phi)))
}
