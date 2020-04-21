#' Rock structure data from an oil well
#' 
#' @description This data comes from lowering a probe into a bore-hole, and taking measurements of the rock structure  as  the  probe  is  lowered. As  the probe moves from one rock strata to another we expect to see an abrupt change in the signal from the measurements.
#' 
#' @source Ruanaidh, Joseph JK O., and William J. Fitzgerald. Numerical Bayesian methods applied to signal processing. Springer Science & Business Media, 2012. \url{https://doi.org/10.1007/978-1-4612-0717-7}
#' @format A numeric vector of 4050 obervations
#' @examples 
#'  
#'  # removing outliers
#'  n = length(oilWell)
#'  h = 32
#'  med = rep(NA, n)
#'  for (i in 1:n) {
#'    index = max(1, i - h):min(n, i + h)
#'    med[i] = median(oilWell[index])
#'  }
#'  residual = (oilWell - med)
#'  
#'  y = oilWell[abs(residual) < 8000]
#'  sigma = sqrt(var(residual[abs(residual) < 8000]))
#'  
#'  # running DeCAFS
#'  res <- DeCAFS(y/sigma)
#'  plot(res, xlab = "time", ylab = "y", type = "l")
#'  abline(v = res$changepoints, col = 4, lty = 3)
"oilWell"

