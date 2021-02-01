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
dataOrnsteinUhlenbech <- function (n = 1e3, y0 = 0, theta = 0, sdEta = 1, sdNu = 1, sdEpsilon = 1, type = c("none", "up", "updown", "rand1"), nbSeg = 20, jumpSize = 10)  {

  f <- scenarioGenerator(n, type = type, nbSeg = nbSeg, jumpSize = jumpSize)
  dw  <- rnorm(n, 0, sdEta)

  nu <- vector(length = n, mode = "double")
  nu[1] <- rnorm(1, 0, sdNu^2/(1-theta^2) )
  for (t in 2:n) {
    nu[t] <-  nu[t-1] - theta * nu[t-1] + sdNu * dw[t-1]
  }

  signal <- f + nu
  y <- f + nu + sdEpsilon

  return(list(y = y, signal = signal, changepoints = which(diff(f) != 0)))
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
                           jumpSize = 5)


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


# generate the F1 dataset
F1df <- mclapply(1:nrow(toSummarize), function(i) {
  p <- toSummarize[i, ]
  print(p)

  fileName <- paste(c("simulations/additional_simulations/resOU/", p, ".RData"), collapse = "")
  if (!file.exists(fileName)) {
    cat("Missing", paste0(Map(paste, names(p), p), collapse = " "), "\n")
    return(NULL)
  } else load(fileName)

  AR1segdfest <- cbind(p$theta,
                       sapply(resar1segEST, function(r)
                         computeF1Score(c(changepoints, N), r$PPSelectedBreaks, 3)) %>% as.numeric, # comptuting the F1
                       sapply(resar1segEST, function(r)
                         computePrecision(c(changepoints, N), r$PPSelectedBreaks, 3)) %>% as.numeric, # comptuting the F1
                       sapply(resar1segEST, function(r)
                         computeRecall(c(changepoints, N), r$PPSelectedBreaks, 3)) %>% as.numeric, # comptuting the F1
                       as.character(p$scenario),
                       "AR1Seg est")
  DeCAFSdfK15 <- cbind(p$theta,
                       sapply(resDeCAFSESTK15, function(r) computeF1Score(c(changepoints, N), c(r$changepoints, N), 3)) %>% as.numeric,
                       sapply(resDeCAFSESTK15, function(r) computePrecision(c(changepoints, N), c(r$changepoints, N), 3)) %>% as.numeric,
                       sapply(resDeCAFSESTK15, function(r) computeRecall(c(changepoints, N), c(r$changepoints, N), 3)) %>% as.numeric,
                       as.character(p$scenario),
                       "DeCAFS est")


  return(rbind(DeCAFSdfK15, AR1segdfest))
}, mc.cores = 6)


F1df <- Reduce(rbind, F1df)
colnames(F1df) <- c("theta", "F1Score", "Precision", "Recall", "Scenario", "Algorithm")
F1df <- as_tibble(F1df) %>% mutate(theta = as.numeric(theta),
                                   F1Score = as.numeric(F1Score),
                                   Precision = as.numeric(Precision),
                                   Recall = as.numeric(Recall))


cbPalette3 <- c("#56B4E9",  "#33cc00", "#434242", "#542354")
scores <- ggplot(F1df,
                 aes(x = theta, y = F1Score, group = Algorithm, by = Algorithm, col = Algorithm)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .05) +
  facet_wrap(~ Scenario) +
  scale_color_manual(values = cbPalette3) +
  xlab(expression(italic(theta)))

scores
ggsave(scores, width = 6, height = 4, units = "in", file = "simulations/outputs/OUF1.pdf", device = "pdf", dpi = "print")


Prec <- ggplot(F1df,
                 aes(x = theta, y = Precision, group = Algorithm, by = Algorithm, col = Algorithm)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .05) +
  facet_wrap(~ Scenario) +
  scale_color_manual(values = cbPalette3) +
  xlab(expression(italic(theta)))

Prec
ggsave(Prec, width = 6, height = 4, units = "in", file = "simulations/outputs/OUPrec.pdf", device = "pdf", dpi = "print")


Rec <- ggplot(F1df,
                 aes(x = theta, y = Recall, group = Algorithm, by = Algorithm, col = Algorithm)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .05) +
  facet_wrap(~ Scenario) +
  scale_color_manual(values = cbPalette3) +
  xlab(expression(italic(theta)))

Rec
ggsave(Rec, width = 6, height = 4, units = "in", file = "simulations/outputs/OURec.pdf", device = "pdf", dpi = "print")
