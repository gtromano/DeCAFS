# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                       #
#   Comparing LAVA oracle vs cross-val to assess which                  #
#   has better perfomances                                              #
#                                                                       #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

library(tidyverse)
library(parallel) # for mclapply
library(Lavash)
library(DeCAFS)

source("simulations/helper_functions.R")

lavaCHANGEPOINT <- function(y, l1penalty, l2penalty) {
  N <- length(y)
  # creating a lower triangular matrix
  L = matrix(1, N, N)
  L[upper.tri(L, diag = T)] <- 0

  K <- 2
  res <- Lavash(L, matrix(y), K, l1penalty, l2penalty, method="profile", Maxiter = 50)
  #return(c(which(res$lava_sparse != 0), N))
  return(list(res = res, est = L %*% res$post_lava, cp = which(res$postlava_sparse != 0)))
}


getLavaPenalty <- function (sdEta, sdNu, N) {
  if (sdEta != 0)
    return( sdNu^2 * (1 / (N * sdEta ^ 2)) )
  else
    return(1e3) # instead of INF
}


REPS <- 96 # number of replicates
N <- 1e3 # lenght of the sequence
CORES <- 16

# range of model parameters
amplitudes <- c(2, 5)
frequencies <- .005
stds <- 1
# jump size
jumpSizes <- c(5)
# scenarios
scenarios <- c("updown")
# segments
nbSegs <- c(4, 20)
# generate a list of simulations
simulations <- expand.grid(amplitude = amplitudes, frequency = frequencies, sd = stds, scenario = scenarios, jumpSize = jumpSizes, nbSeg = nbSegs)


##### FUNCTION FOR RUNNING SIMULATIONS #####

runSim <- function(i, simulations) {
  fileName <- paste(c("simulations/additional_simulations/resLAVAoraclecv/", simulations[i, ], ".RData"), collapse = "")
  if (!file.exists(fileName)) {
    cat("Running ", fileName, "\n")
    #load(fileName)

    p <- simulations[i, ]
    Y <- mclapply(1:REPS, function(r) dataSinusoidal(N, amplitude = p$amplitude, frequency = p$frequency, sd = p$sd, type = as.character(p$scenario), jumpSize = p$jumpSize, nbSeg = p$nbSeg), mc.cores = CORES)

    signal <- lapply(Y, function(r) r$signal)
    y <- lapply(Y, function(r) r$y)
    changepoints <- Y[[1]]$changepoints

    # estimate of the "true" reciprocal of the l2 penalty
    estVariation <- diff(dataSinusoidal(N, amplitude = p$amplitude, frequency = p$frequency, sd = p$sd)$signal) ^ 2 %>% mean %>% sqrt

    # LAVA oracle
    resLAVA <- mclapply(y, lavaCHANGEPOINT, l1penalty = seq(0.01,6,6/50), l2penalty = getLavaPenalty(estVariation, p$sd, N), mc.cores = CORES)

    resLAVACV <- mclapply(y, lavaCHANGEPOINT, l1penalty = seq(0.01,6,6/50), l2penalty = c(0.01, 0.07, 0.2, 0.7, 3,10,60,1000, 2000), mc.cores = CORES)

    save(signal,
         y,
         changepoints,
         resLAVA,
         resLAVACV,
         file = fileName)
  }
}


# running simulations
toSummarize <- simulations %>% filter(nbSeg == 20)


if (T) lapply(1:nrow(toSummarize), runSim, simulations = toSummarize)
