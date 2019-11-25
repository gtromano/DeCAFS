\name{l2FPOP-package}
\alias{l2FPOP-package}
\alias{l2FPOP}
\docType{package}
\title{
  Functional Pruning Optimal Partitioning for ill-posed problems
}
\description{
  Detecting abrupt changes in mean in presence of autocorrelation or random fluctuations.
}
\details{
  Application of the FPOP algorithm with an added l2 penalty for ill-posed problems such as data generated from a Random Walk or an Autoregressive model. The l2fpop function takes now an additional penalty made to cope with the higly correlated data. 
}
\author{
Paul Fearnhead, Guillem Rigaill, Gaetano Romano, Vincent Runge.

Maintainer: Gaetano Romano <g.romano@lancaster.ac.uk>
}
\references{
  To link to the paper
}
\keyword{ package }
\seealso{
  Optional links to other man pages
}
\examples{
  \dontrun{
  set.seed(42)
  Y <- dataRWAR(1000)
  # the l2 parameter is referring to the second penalty
  res <- l2fpop(Y$y, beta = 2 * log(length(Y$y)), lambda = 1, gamma = 1, phi = .98)
  plot(Y$y)
  abline(v = res$changepoints)
  abline(v = Y$changepoints, col = "red", lty = 2, lwd = 2)
  legend("topright", legend = c("estimated", "real"), col = c(1, 2), lty = 1:2)
  }
}
