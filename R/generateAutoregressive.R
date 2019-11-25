.dataAR <- function(n, poisParam = 0.05, meanGap = 10, y0 = 10, gamma = .99, sdY = 1) {
  changepoints <- rpois(n, poisParam)
  cp = which(changepoints > 0)
  mu = cumsum(sample(c(-1, 1), size = n, replace = TRUE) * changepoints * rnorm(n, mean = meanGap))
  y <- .dataAR_c(gamma, y0, mu, rnorm(n, sd = sdY))$z
  return(list(y = y, cp = cp))
}