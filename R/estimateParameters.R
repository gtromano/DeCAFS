require(robust)

estimateParameters = function(y) {
  N   <- length(y)
  yt  <- y[1:(N-2)]
  yt1 <- y[2:(N-1)]
  yt2 <- y[3:N]
  
  cov <- covRob(cbind(yt, yt1, yt2))$cov
  
  phi <- cov[2, 1] / cov[1, 1]
  sdNi <- cov[2, 1] ^ 2 / (-cov[3, 1] + cov[2, 1])
  sdEta <- (cov[1,1] * cov[2,1] ^ 2 - cov[1,1] * cov[3,1] ^ 2 + 2 * cov[2,1] ^ 3)/((-cov[3,1] + cov[2,1])*(cov[2,1] + cov[3,1]))
  return(list(phi = phi, sdEta = sqrt(sdEta), sdNi = sqrt(sdNi)))
}

