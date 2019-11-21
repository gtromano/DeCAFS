estimateParameters = function(y, K = 20) {
  n   <- length(y)

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
  sdEta <- coefficients(model)[1]
  sdNu <- coefficients(model)[2]
  
  # re-estimating the phi
  phi = -(-2 * sdNu + Wk2[1] - sdEta)/(Wk2[1] - sdEta)
  
  return(list(sdEta = as.numeric(sqrt(sdEta)), sdNu = as.numeric(sqrt(sdNu)), phi = as.numeric(phi)))
}
