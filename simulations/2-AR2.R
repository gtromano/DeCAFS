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
dataARp <- function(n = 1e3, poisParam = 0.01, meanGap = 10, phivec = c(.9, -.5), sdNu = 1) {
  changepoints <- rpois(n, poisParam)
  f = cumsum(sample(c(-1, 1), size = n, replace = TRUE) * changepoints * rnorm(n, mean = meanGap))
  mu = f
  
  epsilon <- arima.sim(n = n, list(ar = phivec), sd = sdNu)
  
  y = epsilon + mu
  return(list(y = y, signal = mu, changepoints = which(changepoints > 0)))
}


scenarioGenerator <- function(N, type = c("none", "up", "updown", "rand1"), jumpSize) {
  
  set.seed(42)
  rand1CP <- rpois(20, lambda = 10)
  rand1CP <- rand1CP / max(rand1CP)
  rand1CP <- rand1CP / sum(rand1CP)
  
  set.seed(43)
  rand1Jump <- runif(20, min = -1, max = 1)
  
  type <- match.arg(type)
  switch(
    type,
    none = rep(0, N),
    up = lapply(0:4, function (k)
      rep(k * jumpSize, N * 1 / 5)) %>% unlist,
    updown = lapply(0:4, function(k)
      rep((k %% 2) * jumpSize, N * 1 / 5)) %>% unlist,
    rand1 = c(
      rep(.8 * jumpSize, N * (2 / 5)),
      rep(1.2 * jumpSize, N * (1 / 5)),
      rep(2 * jumpSize, N * (1 / 10)),
      rep(.2 * jumpSize, N * (1 / 10)),
      rep(-.1 * jumpSize, N * (1 / 5))
    )
  )
}


# run the simulation
runSim <- function(i) {
  fileName <- paste(c("simulations/resAR2/", simulations[i, ], ".RData"), collapse = "")
  if (!file.exists(fileName)) {
    cat("Running ", fileName, "\n")
    p <- simulations[i, ]
    
    Z <- scenarioGenerator(N, type = as.character(p$scenario), jumpSize = p$jumpSize)
    
    Y <- lapply(1:REPS, function(r) dataARp(N, poisParam = 0, phivec = c(p$phi, p$phi2),  sd = p$sd))
    signal <- lapply(1:REPS, function(r) Y[[r]]$signal + Z)
    y <- lapply(1:REPS, function(r) Y[[r]]$y + Z)
    changepoints <- which(diff(Z) != 0) + 1
    
    resar1seg <- lapply(y, AR1seg_func, Kmax = 10)
    
    # DeCAFS K 15
    resDeCAFSESTK15 <- lapply(y, DeCAFS)
    
    save(signal, y, changepoints, resar1seg, resDeCAFSESTK15, file = fileName)
  }
}






##### RUNNING SIMULATIONS ####

if (F) mclapply(1:nrow(simulations), runSim, mc.cores = 8)

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
                       sapply(resDeCAFSESTK15, function(r)
                         computeF1Score(c(changepoints, N), c(r$changepoints, N), 3)) %>% as.numeric,
                       as.character(p$scenario), 
                       "DeCAFS est")
  AR1segdf <- cbind(p$phi,
                    p$phi2,
                    p$sd,
                    sapply(resar1seg, function(r)
                      computeF1Score(c(changepoints, N), r$PPSelectedBreaks, 3)) %>% as.numeric,
                    as.character(p$scenario), 
                    "AR1Seg est")
  return(rbind(AR1segdf, DeCAFSdfK15))
}, mc.cores = 4)


F1df <- Reduce(rbind, F1df)

colnames(F1df) <- c("phi", "phi2", "sd", "F1Score",  "Scenario", "Algorithm")
F1df <- as_tibble(F1df) %>% mutate(phi = as.numeric(phi),
                                   phi2 = as.numeric(phi2),
                                   sd = as.numeric(sd),
                                   F1Score = as.numeric(F1Score))

save(F1df, file = "simulations/outputs/F1AR2.RData")
load("simulations/outputs/F1AR2.RData")

cbPalette3 <- c("#56B4E9",  "#33cc00")
scores <- ggplot(F1df %>%  filter(Algorithm != "DeCAFS est (5)" & Algorithm != "DeCAFS est (10)"),
                 aes(x = phi2, y = F1Score, group = Algorithm, color = Algorithm, by = Algorithm)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .05) +
  facet_wrap(. ~ Scenario ) + 
  scale_color_manual(values = cbPalette3) +
  xlab(expression(phi[2]))
scores
ggsave(scores, width = 9, height = 3, units = "in", file = "simulations/outputs/4-missclasAR2.pdf", device = "pdf", dpi = "print")

ggsave(scores, width = 6, height = 4, units = "in", file = "simulations/outputs/4-missclasAR2.pdf", device = "pdf", dpi = "print")
