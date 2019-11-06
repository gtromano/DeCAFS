dataRWAR <- function(n = 1e3, poisParam = 0.01, meanGap = 10, phi = .98, sdEta = 1, sdNi = 1) {
  changepoints <- rpois(n, poisParam)
  f = cumsum(sample(c(-1, 1), size = n, replace = TRUE) * changepoints * rnorm(n, mean = meanGap))
  g = cumsum(rnorm(n, 0, sdEta))
  mu = f + g
  epsilon <- dataAR_c(phi, rnorm(1, 0, sdNi/sqrt(1 - phi^2)), mu, rnorm(n, sd = sdNi))$z
  y = epsilon
  return(list(y = y, f = f, changepoints = which(changepoints > 0)))
}
