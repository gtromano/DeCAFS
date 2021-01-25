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
CORES <- 6

# generate a list of simulations
simulations <- expand.grid(sigmaEta = 1,
                           sigmaNu = 1,
                           theta = seq(0, .99, length.out = 8),
                           scenario = c("none", "up", "updown", "rand1"),
                           jumpSize = 10)


##### FUNCTION FOR RUNNING SIMULATIONS ####
runSim <- function(i, simulations) {
  # here we save the simulation
  fileName <- paste(c("simulations/resRWAR/", simulations[i, ], ".RData"), collapse = "")
  # if file already exist, do not run
  if (!file.exists(fileName)) {
    cat("Running ", fileName, "\n")
    p <- simulations[i, ]

    Y <- mclapply(1:REPS, function(r) dataRWAR(N, phi = p$phi, sdEta = p$sigmaEta, sdNu = p$sigmaNu,  jumpSize = p$jumpSize, type = as.character(p$scenario)), mc.cores = 6)

    signal <- lapply(Y, function(r) r$signal)
    y <- lapply(Y, function(r) r$y)
    changepoints <- Y[[1]]$changepoints


    # this is DeCAFS
    resDeCAFS <- mclapply(y, DeCAFS, beta = (2 * log(N)), modelParam = list(sdEta = p$sigmaEta, sdNu = p$sigmaNu, phi = p$phi), mc.cores = 6)

    # this is fpop
    resfpop <- lapply(y, Fpop, lambda = (2 * (p$sigmaNu^2) * log(N))) # the lambda here is the beta

    # fpop inflated penalty
    if(p$phi != 0)
      resenffpop <- lapply(y, Fpop, lambda = (2 * (p$sigmaNu^2) * (1 + p$phi) / (1 - p$phi) * log(N)))
    else
      resenffpop <- resfpop

    # ar1seg
    resar1seg <- mclapply(y, function(y) AR1seg_func(y, Kmax = 40, rho = p$phi), mc.cores = 6)

    # threshold method
    resThreshold <- mclapply(y, l2Threshold, beta = 2 * log(length(y)), lambda = 1/(p$sigmaEta^2), mc.cores = 6)

    # DeCAFS K 15
    if (p$sigmaEta == 0)
      resDeCAFSESTK15 <- mclapply(y, function(y){
        est <- estimateParameters(y, sdEtaUpper = .0001)
        est$sdEta <- 0
        DeCAFS(y, beta = 2 * log(N), modelParam = est)
      }, mc.cores = 6)
    else
      resDeCAFSESTK15 <- mclapply(y, DeCAFS, mc.cores=6)

    # ar1seg with estimator
    resar1segEST <- mclapply(y, AR1seg_func, Kmax = 40, mc.cores=6)

    # threshold with estimator
    resThresholdEST15 <- mclapply(y, function(y){
      est <- estimateParameters(y)
      l2Threshold(y, beta = 2 * log(N),  lambda =  1/est$sdEta^2)
    }, mc.cores = 6)


    save(
      signal,
      y,
      changepoints,
      resDeCAFS,
      resfpop,
      resenffpop,
      resar1seg,
      resThreshold,
      resDeCAFSESTK15,
      resar1segEST,
      resThresholdEST15,
      file = fileName
    )
  }
}

