# generateRandomWalk = function(n, poisParam = 0.05, avgchange = 5, sd = 1) {
#   changes = rpois(n, poisParam)
#   y = rep(0, n)
#   
#   for (i in 2:n) {
#     y[i] = y[i - 1] + rnorm(1, 0, sd) + sample(c(-1, 1), 1) * changes[i] * rnorm(1, avgchange)
#   }
#   
#   return(list(y = y, cp = which(changes > 0)))
# }

dataRW <- function(n, poisParam = 0.05, meanGap = 10, sdX = 1, sdY = 1) {
  changepoints <- rpois(n, poisParam)
  z <- rnorm(n, 0, sdX) + sample(c(-1,1), size = n, replace = TRUE) * changepoints * rnorm(n, mean = meanGap)
  y <- cumsum(z) + rnorm(n, 0, sdY)
  
  return(list(y = y, cp = which(changepoints > 0)))
}