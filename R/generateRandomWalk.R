.dataRW <- function(n, poisParam = 0.05, meanGap = 10, sdX = 1, sdY = 1) {
  changepoints <- rpois(n, poisParam)
  z <- rnorm(n, 0, sdX) + sample(c(-1,1), size = n, replace = TRUE) * changepoints * rnorm(n, mean = meanGap)
  y <- cumsum(z) + rnorm(n, 0, sdY)
  
  return(list(y = y, cp = which(changepoints > 0)))
}