library(tidyverse)
library(parallel) # for mclapply
library(DeCAFS)
library(Lavash)

source("simulations/helper_functions.R")

lavaCHANGEPOINT <- function(y, l1penalty, l2penalty) {
  N <- length(y)
  # creating a lower triangular matrix
  L = matrix(1, N, N)
  L[upper.tri(L, diag = T)] <- 0
  
  K <- 2
  Lambda1 <- rep(l1penalty, 2)
  Lambda2 <- l2penalty
  res <- Lavash(L, matrix(y), K, Lambda1, Lambda2, method="profile", Maxiter = 50)
  #return(c(which(res$lava_sparse != 0), N))
  return(list(res = res, est = L %*% res$post_lava, cp = which(res$postlava_sparse != 0)))
}


Y <- dataSinusoidal(
  1e3,
  frequency = .005,
  amplitude = 3,
  type = "updown",
  jumpSize = 5,
  nbSeg = 4,
  sd = 1
)


plot(Y$y)

system.time(res <- lavaCHANGEPOINT(Y$y, l1penalty = .1, l2penalty = .01))
lines(res$est, col = 2)
mse(Y$signal, res$est)
computeF1Score(c(res$cp, 1e3), real = c(Y$changepoints, 1e3))


lines(DeCAFS(Y$y)$signal, col = 3)



REPS <- 100 # number of replicates
N <- 1e3 # lenght of the sequence

# range of model parameters
amplitudes <- seq(10, 20, length.out = 5)
frequencies <- seq(10, 100, length.out = 5) / 1e3
stds <- seq(1, 10, length.out = 10)

# jump size
jumpSizes <- c(5)

# scenarios
scenarios <- c("none", "up", "updown", "rand1")

# generate a list of simulations
simulations <- expand.grid(amplitude = amplitudes, frequency = frequencies, sd = stds, scenario = scenarios, jumpSize = jumpSizes)

