# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                           #
#   This first batch of simulations runs simulations for accuracy           #
#   over a range of parameters on the Random Walk AR model                  #
#   The first part of this script creates a framework to run the            #
#   simulations, the second summarises different results, specifically for: #
#    - ranging values of the autocorrelation parameter                      #
#    - different jump sizes                                                 #
#    - different values of the RW noise                                     #
#                                                                           #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

library(tidyverse)
library(parallel) # for mclapply
library(DeCAFS)   
library(AR1seg)
library(fpop)

source("./helper_functions.R")



##### SETTING UP SIMULATIONS #####

# range of model parameters
sigmaEtas <- seq(0, 3, by = .5)
sigmaNus <- seq(1, 5, by = 1)
phis <- seq(0, .99, length.out = 8)

# simulation settings
REPS <- 100 # number of replicates
N = 5e3 # lenght of the sequence

# jump sizes
jumpSizes <- seq(1, 20, by = 3) 

# scenarios (see function scenarioGenerator in helper_functions.R to add more)
scenarios <- c("none", "up", "updown", "rand1")

# generate a list of simulations
simulations <- expand.grid(sigmaEta = sigmaEtas, sigmaNu = sigmaNus, phi = phis, scenario = scenarios, jumpSize = jumpSizes)



##### FUNCTION FOR RUNNING SIMULATIONS ####
runSim <- function(i, simulations) {
  
  # here we save the simulation
  fileName <- paste(c("simulations/resRWAR/", simulations[i, ], ".RData"), collapse = "")
  
  # if file already exist, do not run
  if (!file.exists(fileName)) {
    cat("Running ", fileName, "\n")
    p <- simulations[i, ]
    
    Z <- scenarioGenerator(N, type = as.character(p$scenario), jumpSize = p$jumpSize)
    
    Y <- lapply(1:REPS, function(r) dataRWAR(N, poisParam = 0, phi = p$phi, sdEta = p$sigmaEta, sdNu = p$sigmaNu))
    signal <- lapply(1:REPS, function(r) Y[[r]]$signal + Z)
    y <- lapply(1:REPS, function(r) Y[[r]]$y + Z)
    changepoints <- which(diff(Z) != 0)
    

    # this is DeCAFS
    resDeCAFS <- lapply(y, DeCAFS, beta = (2 * log(N)), modelParam = list(sdEta = p$sigmaEta, sdNu = p$sigmaNu, phi = p$phi))
    
    # this is fpop
    resfpop <- lapply(y, Fpop, lambda = (2 * (p$sigmaNu^2) * log(N))) # the lambda here is the beta
    
    # fpop inflated penalty
    if(p$phi != 0)
      resenffpop <- lapply(y, Fpop, lambda = (2 * (p$sigmaNu^2) * (1 + p$phi) / (1 - p$phi) * log(N)))
    else
      resenffpop <- resfpop
    
    # ar1seg
    resar1seg <- lapply(y, function(y) AR1seg_func(y, Kmax = 40, rho = p$phi))
    
    # threshold method
    resThreshold <- lapply(y, l2Threshold, beta = 2 * log(length(y)), lambda = 1/(p$sigmaEta^2))
    
    # DeCAFS K 15
    if (p$sigmaEta == 0)
      resDeCAFSESTK15 <- resDeCAFSESTK10 <- lapply(y, function(y){
        est <- estimateParameters(y, sdEtaUpper = .0001)
        est$sdEta <- 0
        DeCAFS(y, beta = 2 * log(N), modelParam = est)
      })
    else
      resDeCAFSESTK15 <- lapply(y, DeCAFS)
    
    # ar1seg with estimator
    resar1segEST <- lapply(y, AR1seg_func, Kmax = 40)
    
    # threshold with estimator
    resThresholdEST15 <- lapply(y, function(y){
      est <- estimateParameters(y)
      l2Threshold(y, beta = 2 * log(N),  lambda =  1/est$sdEta^2)
    })
    
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




##### PART 1 - AR ONLY FOR RANGING PHI #####

# selecting the relevant simulations
toSummarize <- simulations %>% filter(sigmaEta == 0, sigmaNu == 2, jumpSize == 10)


# running the simulations (set to T to run)
if (F) mclapply(1:nrow(toSummarize), runSim, simulations = toSummarize, mc.cores = 8)

# generate the F1 dataset
F1df <- mclapply(1:nrow(toSummarize), function(i) {
  p <- toSummarize[i, ]
  print(p)
  
  fileName <- paste(c("simulations/resRWAR/", p, ".RData"), collapse = "")
  if (!file.exists(fileName)) {
    cat("Missing", paste0(Map(paste, names(p), p), collapse = " "), "\n")
    return(NULL)
  } else load(fileName)
  
  DeCAFSdf <- cbind(p$phi,
                    sapply(resDeCAFS, function(r) computeF1Score(c(changepoints, N), c(r$changepoints, N), 3)) %>% as.numeric,
                    as.character(p$scenario), 
                    "DeCAFS")
  fpopdf <- cbind(p$phi,
                  sapply(resfpop, function(r) computeF1Score(c(changepoints, N), r$t.est, 3)) %>% as.numeric, 
                  as.character(p$scenario), 
                  "fpop")
  enffpopdf <- cbind(p$phi,
                     sapply(resenffpop, function(r) computeF1Score(c(changepoints, N), r$t.est, 3)) %>% as.numeric,
                     as.character(p$scenario), 
                     "fpop Inf")
  AR1segdf <- cbind(p$phi,
                    sapply(resar1seg, function(r)
                      computeF1Score(c(changepoints, N), r$PPSelectedBreaks, 3)) %>% as.numeric, # comptuting the F1
                    as.character(p$scenario), 
                    "AR1Seg")
  
  AR1segdfest <- cbind(p$phi,
                       sapply(resar1segEST, function(r)
                         computeF1Score(c(changepoints, N), r$PPSelectedBreaks, 3)) %>% as.numeric, # comptuting the F1
                       as.character(p$scenario), 
                       "AR1Seg est")
  DeCAFSdfK15 <- cbind(p$phi,
                       sapply(resDeCAFSESTK15, function(r) computeF1Score(c(changepoints, N), c(r$changepoints, N), 3)) %>% as.numeric,
                       as.character(p$scenario), 
                       "DeCAFS est")
  
  return(rbind(DeCAFSdf, fpopdf, enffpopdf, AR1segdf, DeCAFSdfK15, AR1segdfest))
}, mc.cores = 4)


F1df <- Reduce(rbind, F1df)
colnames(F1df) <- c("phi", "F1Score", "Scenario", "Algorithm")
F1df <- as_tibble(F1df) %>% mutate(phi = as.numeric(phi),
                                   F1Score = as.numeric(F1Score))

# save dataset
save(F1df, file = "simulations/outputs/F1AR.RData")

# load dataset
load("simulations/outputs/F1AR.RData")

cbPalette <- c("#0072B2", "#56B4E9", "#009E73", "#33cc00", "#E69F00", "#CC79A7")
scores <- ggplot(F1df %>% filter(Algorithm != "DeCAFS est (5)" & Algorithm != "DeCAFS est (10)"),
                 aes(x = phi, y = F1Score, group = Algorithm, by = Algorithm, col = Algorithm)) +
  geom_vline(xintercept = unique(simulations$phi)[7], col = "grey", lty = 2) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .05) +
  facet_wrap(~ Scenario) +
  scale_color_manual(values = cbPalette) + 
  xlab(expression(italic(phi)))

scores

ggsave(scores,  width = 7, height = 7, file = "simulations/outputs/1-AR.pdf", device = "pdf", dpi = "print")




##### PART 2 - AR ONLY FOR VARIOUS JUMP SIZES #####

# selecting the relevant simulations
toSummarize <- simulations %>% filter(sigmaEta == 0, sigmaNu == 2, phi == unique(simulations$phi)[7])


# running the simulations
if (F) mclapply(1:nrow(toSummarize), runSim, simulations = toSummarize, SIMID = "", mc.cores = 8)


F1df <- mclapply(1:nrow(toSummarize), function(i) {
  p <- toSummarize[i, ]
  print(p)
  
  fileName <- paste(c("simulations/resRWAR/", p, ".RData"), collapse = "")
  if (!file.exists(fileName)) {
    cat("Missing", paste0(Map(paste, names(p), p), collapse = " "), "\n")
    return(NULL)
  } else load(fileName)
  
  DeCAFSdf <- cbind(p$jumpSize,
                    sapply(resDeCAFS, function(r) computeF1Score(c(changepoints, N), c(r$changepoints, N), 3)) %>% as.numeric,
                    as.character(p$scenario), 
                    "DeCAFS")
  fpopdf <- cbind(p$jumpSize,
                  sapply(resfpop, function(r) computeF1Score(c(changepoints, N), r$t.est, 3)) %>% as.numeric, 
                  as.character(p$scenario), 
                  "fpop")
  enffpopdf <- cbind(p$jumpSize,
                     sapply(resenffpop, function(r) computeF1Score(c(changepoints, N), r$t.est, 3)) %>% as.numeric,
                     as.character(p$scenario), 
                     "fpop Inf")
  AR1segdf <- cbind(p$jumpSize,
                    sapply(resar1seg, function(r)
                      computeF1Score(c(changepoints, N), r$PPSelectedBreaks, 3)) %>% as.numeric, # comptuting the F1
                    as.character(p$scenario), 
                    "AR1Seg")
  AR1segdfest <- cbind(p$jumpSize,
                       sapply(resar1segEST, function(r)
                         computeF1Score(c(changepoints, N), r$PPSelectedBreaks, 3)) %>% as.numeric, # comptuting the F1
                       as.character(p$scenario), 
                       "AR1Seg est")
  DeCAFSdfK15 <- cbind(p$jumpSize,
                       sapply(resDeCAFSESTK15, function(r) computeF1Score(c(changepoints, N), c(r$changepoints, N), 3)) %>% as.numeric,
                       as.character(p$scenario), 
                       "DeCAFS est")
  
  return(rbind(DeCAFSdf, fpopdf, enffpopdf, AR1segdf,  DeCAFSdfK15, AR1segdfest))
}, mc.cores = 4)


F1df <- Reduce(rbind, F1df)

colnames(F1df) <- c("jumpSize", "F1Score", "Scenario", "Algorithm")
F1df <- as_tibble(F1df) %>% mutate(jumpSize = as.numeric(jumpSize),
                                   F1Score = as.numeric(F1Score))

save(F1df, file = "simulations/outputs/F1ARJumpSize.RData")

# load dataset
load("simulations/outputs/F1ARJumpSize.RData")

cbPalette <- c("#0072B2", "#56B4E9", "#009E73", "#33cc00", "#E69F00", "#CC79A7")
scoresJump <- ggplot(F1df %>% filter(Algorithm != "DeCAFS est (5)" & Algorithm != "DeCAFS est (10)"),
                     aes(x = jumpSize, y = F1Score, group = Algorithm, by = Algorithm, col = Algorithm)) +
  geom_vline(xintercept = 10, col = "grey", lty = 2) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .5) +
  facet_wrap(~ Scenario) + 
  scale_color_manual(values = cbPalette) + 
  xlab("Jump Size")
scoresJump



##### PART 3 - RWAR FOR RANGING SDETA #####


# selecting the simulations with sd_nu = 2
toSummarize <- simulations %>% filter(phi == (simulations$phi %>% unique())[7], sigmaNu == 2, jumpSize == 10)

if (F) mclapply(1:nrow(toSummarize), runSim, simulations = toSummarize, SIMID = "", mc.cores = 8)


F1df <- mclapply(1:nrow(toSummarize), function(i) {
  p <- toSummarize[i, ]
  print(p)
  
  fileName <- paste(c("simulations/resRWAR/", p, ".RData"), collapse = "")
  if (!file.exists(fileName)) {
    cat("Missing", paste0(Map(paste, names(p), p), collapse = " "), "\n")
    return(NULL)
  } else load(fileName)
  
  DeCAFSdf <- cbind(p$sigmaEta,
                    sapply(resDeCAFS, function(r) computeF1Score(c(changepoints, N), c(r$changepoints, N), 3)) %>% as.numeric,
                    as.character(p$scenario), 
                    "DeCAFS")
  fpopdf <- cbind(p$sigmaEta,
                  sapply(resfpop, function(r) computeF1Score(c(changepoints, N), r$t.est, 3)) %>% as.numeric, 
                  as.character(p$scenario), 
                  "fpop")
  enffpopdf <- cbind(p$sigmaEta,
                     sapply(resenffpop, function(r) computeF1Score(c(changepoints, N), r$t.est, 3)) %>% as.numeric,
                     as.character(p$scenario), 
                     "fpop Enf")
  AR1segdf <- cbind(p$sigmaEta,
                    sapply(resar1seg, function(r)
                      computeF1Score(c(changepoints, N), r$PPSelectedBreaks, 3)) %>% as.numeric, # comptuting the F1
                    as.character(p$scenario), 
                    "AR1Seg")
  AR1segdfest <- cbind(p$sigmaEta,
                       sapply(resar1segEST, function(r)
                         computeF1Score(c(changepoints, N), r$PPSelectedBreaks, 3)) %>% as.numeric, # comptuting the F1
                       as.character(p$scenario), 
                       "AR1Seg est")
  DeCAFSdfK15 <- cbind(p$sigmaEta,
                       sapply(resDeCAFSESTK15, function(r) computeF1Score(c(changepoints, N), c(r$changepoints, N), 3)) %>% as.numeric,
                       as.character(p$scenario), 
                       "DeCAFS est")
  
  return(rbind(DeCAFSdf, fpopdf, enffpopdf, AR1segdf, DeCAFSdfK15, AR1segdfest))
}, mc.cores = 4)



F1df <- Reduce(rbind, F1df)


colnames(F1df) <- c("SigmaEta", "F1Score", "Scenario", "Algorithm")
F1df <- as_tibble(F1df) %>% mutate(SigmaEta = as.numeric(SigmaEta),
                                   F1Score = as.numeric(F1Score))


save(F1df, file = "simulations/outputs/F1RWAR.RData")

load("simulations/outputs/F1RWAR.RData")

cbPaletteEDIT <- c("#56B4E9", "#009E73", "#33cc00", "#E69F00", "#CC79A7")

scoresRWAR <- ggplot(F1df %>% filter(Algorithm != "DeCAFS est (5)" & Algorithm != "DeCAFS est (10)" & Algorithm != "AR1Seg"),
                     aes(x = SigmaEta, y = F1Score, group = Algorithm, by = Algorithm, col = Algorithm)) +
  geom_vline(xintercept = 0, col = "grey", lty = 2) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.08) +
  facet_wrap(~ Scenario) + 
  scale_color_manual(values = cbPaletteEDIT) + 
  xlab(expression(sigma[eta]))

scoresRWAR



#### COMPLETE PLOT ####

library(ggpubr)

# getting legend out
fullLegend <- get_legend(scores)

outPlot <-
  ggarrange(
    scores,
    scoresJump,
    scoresRWAR,
    as_ggplot(fullLegend),
    labels = c("A", "B", "C", ""),
    legend = "none"
  )
ggsave(outPlot, width = 10, height = 10, file = "simulations/outputs/1-RWARcomplete.pdf", device = "pdf", dpi = "print")

