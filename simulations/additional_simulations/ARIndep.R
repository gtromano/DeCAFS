# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                        #
#   AR w\ independence across various segments           #
#                                                        #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


library(tidyverse)
library(parallel) # for mclapply
library(DeCAFS)
library(AR1seg)

source("simulations/helper_functions.R")


dataARIndep <- function (n = 1e3, phi = 0, sdNu = 1, type = c("none", "up", "updown", "rand1"), nbSeg = 20, jumpSize = 10)  {

  f <- scenarioGenerator(n, type = type, nbSeg = nbSeg, jumpSize = jumpSize)
  changepoints <- which(diff(f) != 0)

  out <- NULL

  for (seglen in diff(c(0, changepoints, n))) {
    y <- dataRWAR(n = seglen, phi = phi, sdNu = 2, sdEta = 0)$y
    out <- c(out, y)
  }

  out <- out + f

  return(list(y = out, signal = f, changepoints = changepoints))
}



# simulation settings
REPS <- 100 # number of replicates
N <- 5e3 # lenght of the sequence
CORES <- 6

# generate a list of simulations
simulations <- expand.grid(sigmaNu = 2,
                           phi = seq(0, .95, length.out = 8),
                           scenario = c("up", "updown", "rand1"),
                           jumpSize = 10)


##### FUNCTION FOR RUNNING SIMULATIONS ####
runSim <- function(i, simulations) {
  # here we save the simulation
  fileName <- paste(c("simulations/additional_simulations/resARIndep/", simulations[i, ], ".RData"), collapse = "")
  # if file already exist, do not run
  if (!file.exists(fileName)) {
    cat("Running ", fileName, "\n")
    p <- simulations[i, ]

    Y <- mclapply(1:REPS, function(r) dataARIndep(N, phi = p$phi, sdNu = p$sigmaNu, jumpSize = p$jumpSize, type = as.character(p$scenario)), mc.cores = CORES)

    signal <- lapply(Y, function(r) r$signal)
    y <- lapply(Y, function(r) r$y)
    changepoints <- Y[[1]]$changepoints

    #DeCAFS K 15
    resDeCAFSESTK15 <- mclapply(y, DeCAFS, mc.cores = CORES)

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


# selecting the relevant simulations
toSummarize <- simulations

# running the simulations (set to T to run)
if (T) lapply(1:nrow(toSummarize), runSim, simulations = toSummarize)




# generate the F1 dataset
F1df <- mclapply(1:nrow(toSummarize), function(i) {
  p <- toSummarize[i, ]
  print(p)

  fileName <- paste(c("simulations/additional_simulations/resARIndep/", p, ".RData"), collapse = "")
  if (!file.exists(fileName)) {
    cat("Missing", paste0(Map(paste, names(p), p), collapse = " "), "\n")
    return(NULL)
  } else load(fileName)

  AR1segdfest <- cbind(p$phi,
                       sapply(resar1segEST, function(r)
                         computeF1Score(c(changepoints, N), r$PPSelectedBreaks, 3)) %>% as.numeric, # comptuting the F1
                       sapply(resar1segEST, function(r)
                         computePrecision(c(changepoints, N), r$PPSelectedBreaks, 3)) %>% as.numeric, # comptuting the F1
                       sapply(resar1segEST, function(r)
                         computeRecall(c(changepoints, N), r$PPSelectedBreaks, 3)) %>% as.numeric, # comptuting the F1
                       as.character(p$scenario),
                       "AR1Seg est")
  DeCAFSdfK15 <- cbind(p$phi,
                       sapply(resDeCAFSESTK15, function(r) computeF1Score(c(changepoints, N), c(r$changepoints, N), 3)) %>% as.numeric,
                       sapply(resDeCAFSESTK15, function(r) computePrecision(c(changepoints, N), c(r$changepoints, N), 3)) %>% as.numeric,
                       sapply(resDeCAFSESTK15, function(r) computeRecall(c(changepoints, N), c(r$changepoints, N), 3)) %>% as.numeric,
                       as.character(p$scenario),
                       "DeCAFS est")


  return(rbind(DeCAFSdfK15, AR1segdfest))
}, mc.cores = 6)


F1df <- Reduce(rbind, F1df)
colnames(F1df) <- c("phi", "F1Score", "Precision", "Recall", "Scenario", "Algorithm")
F1df <- as_tibble(F1df) %>% mutate(phi = as.numeric(phi),
                                   F1Score = as.numeric(F1Score),
                                   Precision = as.numeric(Precision),
                                   Recall = as.numeric(Recall))


cbPalette3 <- c("#56B4E9",  "#33cc00", "#434242", "#542354")
scores <- ggplot(F1df,
                 aes(x = phi, y = F1Score, group = Algorithm, by = Algorithm, col = Algorithm)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .05) +
  facet_wrap(~ Scenario) +
  scale_color_manual(values = cbPalette3) +
  xlab(expression(italic(phi)))

scores
ggsave(scores, width = 6, height = 4, units = "in", file = "simulations/outputs/ARIndepF1.pdf", device = "pdf", dpi = "print")


scores <- ggplot(F1df,
                 aes(x = phi, y = Precision, group = Algorithm, by = Algorithm, col = Algorithm)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .05) +
  facet_wrap(~ Scenario) +
  scale_color_manual(values = cbPalette3) +
  xlab(expression(italic(phi)))

ggsave(scores, width = 6, height = 4, units = "in", file = "simulations/outputs/ARIndepPrec.pdf", device = "pdf", dpi = "print")

scores <- ggplot(F1df,
                 aes(x = phi, y = Recall, group = Algorithm, by = Algorithm, col = Algorithm)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .05) +
  facet_wrap(~ Scenario) +
  scale_color_manual(values = cbPalette3) +
  xlab(expression(italic(phi)))

ggsave(scores, width = 6, height = 4, units = "in", file = "simulations/outputs/ARIndepRecall.pdf", device = "pdf", dpi = "print")
