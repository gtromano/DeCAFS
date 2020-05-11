#' DeCAFS Plotting
#' 
#' DeCAFS output plotting method. 
#'
#' @param x the output object from a DeCAFS call
#' @param ... Additional graphical parameters to be passed down to the plot function
#'
#' @return An R plot
#' 
#' @export 
#' plot.DeCAFSout <- function(object,...) { }
#' @examples
#' set.seed(42)
#' Y = dataRWAR(n = 1e3, poisParam = .01, meanGap = 15, phi = .5, sdEta = 1, sdNu = 3)
#' res = DeCAFS(Y$y)
#' plot(res, type = "l")
#' 


plot.DeCAFSout <- function(x, ...) {
  data <- x$data
  plot(data, ...)
  
  cps <- c(0, x$changepoints, length(x$data)) 
  for(i in 1:(length(cps) - 1))
    lines((cps[i] + 1):cps[i + 1], x$signal[(cps[i] + 1):cps[i + 1]], col = 2, lwd = 2)
  
}
