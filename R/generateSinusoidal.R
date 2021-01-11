#' Generating data from a sinusoidal model with changes
#' 
#' This function generates a sequence of observation from a sinusoidal model with changes. This can be used as an example 
#' for model misspecification. 
#' 
#' @param n The length of the sequence of observations.
#' @param poisParam A poisson parameter regulating the probability of observing a change.
#' @param meanGap The average magnitude of a change.
#' @param amplitude The amplitude of the sinusoid
#' @param frequency The angular frequency of the sinusoid
#' @param phase where the signal starts at time t = 0
#' @param sd standard deviation of the noise added on top of the signal
#'
#' @return A list containing:
#' \describe{
#' \item{\code{y}}{the data sequence,}
#' \item{\code{signal}}{the underlying signal without the noise,}
#' \item{\code{changepoints}}{the changepoint locations}
#' }
#'
#' @export
#'
#' @examples
#' Y <- dataSinusoidal(
#'   1e4,
#'   frequency = 1 / 1e3,
#'   amplitude = 10,
#'   type = "updown",
#'   jumpSize = 4,
#'   nbSeg = 4
#' )
#' res <- DeCAFS(Y$y)
#' plot(res, col = "grey")
#' lines(Y$signal, col = "blue", lwd = 2, lty = 2)
#' abline(v = res$changepoints, col = 2)
#' abline(v = Y$changepoints, col = 4, lty = 2)


dataSinusoidal <- function(n, poisParam = 0.01, amplitude = 1, frequency = 0.001, phase = 0, sd = 1,  type = c("none", "up", "updown", "rand1"), nbSeg = 20, jumpSize = 1) {
  changepoints <- rpois(n, poisParam)
  f <- scenarioGenerator(n, type = type, nbSeg = nbSeg, jumpSize = jumpSize)
  g <- amplitude * sin(2 * pi * frequency * 1:n + phase)
  mu <- f + g
  y <- mu + rnorm(n, sd = sd)
  return(list(y = y, signal = mu, changepoints =  which(diff(f) != 0)))
}
