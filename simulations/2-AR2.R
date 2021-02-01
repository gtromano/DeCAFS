# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                       #
#   This second batch of simulations runs simulations for accuracy      #
#   over a range of parameters on the misspecified process AR2          #
#   The first part of this script creates a framework to run the        #
#   simulations, the second summarises different results                #
#                                                                       #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


library(tidyverse)
library(parallel) # for mclapply
library(DeCAFS)   
library(AR1seg)
library(fpop)

source("simulations/helper_functions.R")


##### SETTING UP SIMULATIONS #####

# simulation settings
REPS <- 100 # number of replicates
N = 5e3 # lenght of the sequence
CORES = 6

# range of model parameters
phi <- .3
phi2 <- seq(-.01, -.99, length.out = 5)
stds <- 2

# jump size
jumpSizes <- c(20)


# scenarios
scenarios <- c("up", "updown", "rand1", "none")

# generate a list of simulations
simulations <- expand.grid(phi =phi, phi2 = phi2, sd = stds, scenario = scenarios, jumpSize = jumpSizes)



##### FUNCTIONS FOR RUNNING SIMULATIONS ####

# generate data from AR2
dataARp <- function(n = 1e3, phivec = c(.9, -.5), sdNu = 1, type = c("none", "up", "updown", "rand1"), nbSeg = 20, jumpSize = 1) {
  mu <- scenarioGenerator(n, type = type, nbSeg = nbSeg, jumpSize = jumpSize)
  epsilon <- arima.sim(n = n, list(ar = phivec), sd = sdNu)
  y <- epsilon + mu
  return(list(y = y, signal = mu, changepoints =  which(diff(mu) != 0)))
}


# run the simulation
runSim <- function(i) {
  fileName <- paste(c("simulations/resAR2/", simulations[i, ], ".RData"), collapse = "")
  if (!file.exists(fileName)) {
    cat("Running ", fileName, "\n")
    p <- simulations[i, ]


    Y <- mclapply(1:REPS, function(r) dataARp(N, phivec = c(p$phi, p$phi2),  sdNu = p$sd, jumpSize = p$jumpSize, type = as.character(p$scenario)), mc.cores = CORES)
    signal <- lapply(Y, function(r) r$signal)
    y <- lapply(Y, function(r) r$y)
    changepoints <- Y[[1]]$changepoints


    #DeCAFS K 15
    resDeCAFSESTK15 <- mclapply(y, DeCAFS, mc.cores = CORES)

    # ar1seg with estimator
    resar1seg <- mclapply(y, AR1seg_func, Kmax = 40, mc.cores=CORES)



    save(signal, y, changepoints, resar1seg, resDeCAFSESTK15, file = fileName)
  }
}



##### RUNNING SIMULATIONS ####

if (T) lapply(1:nrow(simulations), runSim)

toSummarize <- simulations

# summarizing
F1df <- mclapply(1:nrow(toSummarize), function(i) {
  p <- toSummarize[i, ]
  print(p)
  
  fileName <- paste(c("simulations/resAR2/", p, ".RData"), collapse = "")
  if (!file.exists(fileName)) {
    cat("Missing", paste0(Map(paste, names(p), p), collapse = " "), "\n")
    return(NULL)
  } else load(fileName)
  
  DeCAFSdfK15 <- cbind(p$phi,
                       p$phi2,
                       p$sd,
                       sapply(resDeCAFSESTK15, function(r) computeF1Score(c(changepoints, N), c(r$changepoints, N), 3)) %>% as.numeric,
                       sapply(resDeCAFSESTK15, function(r) computePrecision(c(changepoints, N), c(r$changepoints, N), 3)) %>% as.numeric,
                       sapply(resDeCAFSESTK15, function(r) computeRecall(c(changepoints, N), c(r$changepoints, N), 3)) %>% as.numeric,
                       as.character(p$scenario), 
                       "DeCAFS est")
  AR1segdf <- cbind(p$phi,
                    p$phi2,
                    p$sd,
                    sapply(resar1seg, function(r)
                      computeF1Score(c(changepoints, N), r$PPSelectedBreaks, 3)) %>% as.numeric,
                     sapply(resar1seg, function(r)
                       computePrecision(c(changepoints, N), r$PPSelectedBreaks, 3)) %>% as.numeric, # comptuting the F1
                     sapply(resar1seg, function(r)
                       computeRecall(c(changepoints, N), r$PPSelectedBreaks, 3)) %>% as.numeric, # comptuting the F1
                    as.character(p$scenario),
                    "AR1Seg est")
  return(rbind(AR1segdf, DeCAFSdfK15))
}, mc.cores = CORES)


F1df <- Reduce(rbind, F1df)

colnames(F1df) <- c("phi", "phi2", "sd", "F1Score", "Precision",  "Recall", "Scenario", "Algorithm")
F1df <- as_tibble(F1df) %>% mutate(phi = as.numeric(phi),
                                   phi2 = as.numeric(phi2),
                                   sd = as.numeric(sd),
                                   F1Score = as.numeric(F1Score),
                                   Precision = as.numeric(Precision),
                                   Recall = as.numeric(Recall))

save(F1df, file = "simulations/outputs/F1AR2.RData")
load("simulations/outputs/F1AR2.RData")

cbPalette3 <- c("#56B4E9",  "#33cc00")
scores <- ggplot(F1df,
                 aes(x = phi2, y = F1Score, group = Algorithm, color = Algorithm, by = Algorithm)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .05) +
  facet_wrap(. ~ Scenario ) + 
  scale_color_manual(values = cbPalette3) +
  xlab(expression(phi[2]))
scores
ggsave(scores, width = 6, height = 4, units = "in", file = "simulations/outputs/4-missclasAR2.pdf", device = "pdf", dpi = "print")


Prec <- ggplot(F1df,
                 aes(x = phi2, y = Precision, group = Algorithm, color = Algorithm, by = Algorithm)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .05) +
  facet_wrap(. ~ Scenario ) +
  scale_color_manual(values = cbPalette3) +
  xlab(expression(phi[2]))
Prec
ggsave(Prec, width = 6, height = 4, units = "in", file = "simulations/outputs/4-missclasAR2Prec.pdf", device = "pdf", dpi = "print")


Rec <- ggplot(F1df,
                 aes(x = phi2, y = Recall, group = Algorithm, color = Algorithm, by = Algorithm)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .05) +
  facet_wrap(. ~ Scenario ) +
  scale_color_manual(values = cbPalette3) +
  xlab(expression(phi[2]))
Rec
ggsave(Rec, width = 6, height = 4, units = "in", file = "simulations/outputs/4-missclasAR2Rec.pdf", device = "pdf", dpi = "print")
