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
#' @return The sequence of observation, the sinusoid signal at the basis, the changepoints location. 
#' @export
#'
#' @examples
#' Y <- dataSinusoidal(1e4, poisParam = .0005, meanGap = 5, frequency = 2 * pi / 1e3, amplitude = 10, sd = 2)
#' res <- l2fpop(Y$y)
#' par(mfrow = c(2, 1))
#' plot(Y$y, col = "grey")
#' lines(res$signal, col = "blue", lwd = 2)
#' lines(Y$signal, col = "red", lwd = 2, lty = 2)
#' abline(v = res$changepoints, col = 4)
#' abline(v = Y$changepoints, col = 2, lty = 2)
#' plot(Y$y[-1] - Y$y[-1e4])
#' abline(v = Y$changepoints, col = 2, lty = 2)
#' par(mfrow = c(1, 1))

dataSinusoidal <- function(n, poisParam = 0.01, meanGap = 10, amplitude = 1, frequency = 1, phase = 0, sd = 1) {
  changepoints <- rpois(n, poisParam)
  f <- cumsum(sample(c(-1, 1), size = n, replace = TRUE) * changepoints * rnorm(n, mean = meanGap))
  g <- amplitude * sin(frequency * 1:n + phase)
  mu <- f + g
  y <- mu + rnorm(n, sd = sd)
  return(list(y = y, signal = mu, changepoints = which(changepoints > 0)))
}
