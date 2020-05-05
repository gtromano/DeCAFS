#' Generate a Random Walk + AR realization
#'
#' Generate a Realization from the RWAR model (check the references for further details).
#' \deqn{y_t = \mu_t + \epsilon_t}
#' where 
#' \deqn{\mu_t = \mu_{t-1} + \eta_t + \delta_t, \quad \eta_t \sim N(0, \sigma_\eta^2), \ \delta_t \ \in R}
#' and
#' \deqn{\epsilon_t = \phi \epsilon_{t-1} + \nu_t \quad \nu_t \sim N(0, \sigma_\nu^2)}
#' 
#' 
#' @param n The length of the sequence of observations.
#' @param poisParam A poisson parameter regulating the probability of observing a change.
#' @param meanGap The average magnitude of a change.
#' @param phi The autocorrelation parameter \eqn{\phi}
#' @param sdEta The standard deviation of the Random Walk Component on the signal drift
#' @param sdNu The standard deviation of the Autocorrelated noise
#'
#' @return A list containing:
#' \describe{
#' \item{\code{y}}{the data sequence,}
#' \item{\code{signal}}{the underlying signal without the superimposed AR(1) noise,}
#' \item{\code{changepoints}}{the changepoint locations}
#' }
#' 
#' 
#' @export
#' 
#' @references Romano, G., Rigaill, G., Runge, V., Fearnhead, P. Detecting Abrupt Changes in the Presence of Local Fluctuations and Autocorrelated Noise. arXiv preprint \url{https://arxiv.org/abs/2005.01379} (2020).
#'
#' @examples
#' library(ggplot2)
#' set.seed(42)
#' Y = dataRWAR(n = 1e3, poisParam = .01, meanGap = 15, phi = .5, sdEta = 3, sdNu = 1)
#' y = Y$y
#' ggplot(data.frame(t = 1:length(y), y), aes(x = t, y = y)) +
#'   geom_point() +
#'   geom_vline(xintercept = Y$changepoints, col = 4,  lty = 3)


dataRWAR <- function(n = 1e3, poisParam = 0.01, meanGap = 10, phi = 0, sdEta = 0, sdNu = 1) {
  changepoints <- rpois(n, poisParam)
  f = cumsum(sample(c(-1, 1), size = n, replace = TRUE) * changepoints * rnorm(n, mean = meanGap))
  g = cumsum(rnorm(n, 0, sdEta))
  mu = f + g
  epsilon <- .dataAR_c(phi, rnorm(1, 0, sdNu/sqrt(1 - phi^2)), mu, rnorm(n, sd = sdNu))$z
  y = epsilon
  return(list(y = y, signal = mu, changepoints = which(changepoints > 0) - 1))
}
