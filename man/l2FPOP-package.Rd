\name{l2FPOP-package}
\alias{l2FPOP-package}
\alias{l2FPOP}
\docType{package}
\title{
  Functional Pruning Optimal Partitioning for ill-posed problems
}
\description{
  A more detailed description of what the package does. A length
  of about one to five lines is recommended.
}
\details{
  Application of the FPOP algorithm with an added l2 penalty for ill-posed problems such as data generated from a Random Walk or an Autoregressive model. The l2fpop function takes now an additional penalty tailored to cope with the higly correlated data. 
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
  Y = generateRandomWalk(1000, poisParam = 0.005, avgchange = 20, sd = 1) # generating a random walk
  y = Y$y + rnorm(length(Y$y), sd = 1) # adding some noise
  # the l2 parameter is referring to the second penalty
  res = fpopl2(y, l0penalty = 4 * log(length(y)), l2penalty = 1)
  }
}
