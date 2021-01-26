# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                       #
#   A minor simulation on the Ornstein - Uhlenbeck process              #
#   with changes                                                        #
#                                                                       #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

library(tidyverse)
library(parallel) # for mclapply
library(DeCAFS)
library(AR1seg) # not available for R 4.x

source("simulations/helper_functions.R")

# data generating function
dataOrnsteinUhlenbech <- function (n = 1e3, y0 = 0, theta = 0, sdEta = 1, sdNu = 1, type = c("none", "up", "updown", "rand1"), nbSeg = 20, jumpSize = 10)  {

  f <- scenarioGenerator(n, type = type, nbSeg = nbSeg, jumpSize = jumpSize)
  dw  <- rnorm(n, 0, sdEta)
  out <- vector(length = n, mode = "double")
  out[1] <- y0
  for (t in 2:n) {
    out[t] = out[t-1] + theta * (f[t-1] - out[t-1]) + sdNu * dw[t-1]
  }

  return(list(y = out, signal = f, changepoints = which(diff(f) != 0)))
}

# simulation settings
REPS <- 100 # number of replicates
N <- 5e3 # lenght of the sequence
CORES <- 16

# generate a list of simulations
simulations <- expand.grid(sigmaEta = 1,
                           sigmaNu = 1,
                           theta = seq(0, .99, length.out = 8),
                           scenario = c("none", "up", "updown", "rand1"),
                           jumpSize = 10)


##### FUNCTION FOR RUNNING SIMULATIONS ####
runSim <- function(i, simulations) {
  # here we save the simulation
  fileName <- paste(c("simulations/additional_simulations/resOU/", simulations[i, ], ".RData"), collapse = "")
  # if file already exist, do not run
  if (!file.exists(fileName)) {
    cat("Running ", fileName, "\n")
    p <- simulations[i, ]

    Y <- mclapply(1:REPS, function(r) dataOrnsteinUhlenbech(N, theta = p$theta, sdEta = p$sigmaEta, sdNu = p$sigmaNu,  jumpSize = p$jumpSize, type = as.character(p$scenario)), mc.cores = CORES)

    signal <- lapply(Y, function(r) r$signal)
    y <- lapply(Y, function(r) r$y)
    changepoints <- Y[[1]]$changepoints

    # DeCAFS K 15
    if (p$sigmaEta == 0)
      resDeCAFSESTK15 <- mclapply(y, function(y){
        est <- estimateParameters(y, sdEtaUpper = .0001)
        est$sdEta <- 0
        DeCAFS(y, beta = 2 * log(N), modelParam = est)
      }, mc.cores = CORES)
    else
      resDeCAFSESTK15 <- mclapply(y, DeCAFS, mc.cores=CORES)

    # ar1seg with estimator
    resar1segEST <- mclapply(y, AR1seg_func, Kmax = 40, mc.cores=CORES)

    save(
      signal,
      y,
      changepoints,
      resDeCAFSESTK15,
      resar1segEST,
      file = fileName
    )
  }
}




##### RUNNING SIMULATIONS FOR RANGING THETA #####

# selecting the relevant simulations
toSummarize <- simulations

# running the simulations (set to T to run)
if (T) lapply(1:nrow(toSummarize), runSim, simulations = toSummarize)
