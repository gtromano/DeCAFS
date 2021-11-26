# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                       #
#   This second batch of simulations runs simulations for accuracy      #
#   over a range of parameters on the misspecified Sinusoidal process   #
#   The first part of this script creates a framework to run the        #
#   simulations, the second summarises different results                #
#                                                                       #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


library(tidyverse)
library(parallel) # for mclapply
library(DeCAFS)   
library(AR1seg)
library(fpop)
library(ggpubr)

source("simulations/helper_functions.R")


##### SETTING UP SIMULATIONS #####

REPS <- 100 # number of replicates
N = 5e3 # lenght of the sequence
CORES = 6


# range of model parameters
amplitudes <- seq(10, 20, length.out = 5)
frequencies <- seq(10, 100, length.out = 5) / (1e4 * 2 * pi)
stds <- seq(1, 10, length.out = 10)

# jump size
jumpSizes <- c(5)

# scenarios
scenarios <- c("none", "up", "updown", "rand1")

# generate a list of simulations
simulations <- expand.grid(amplitude = amplitudes, frequency = frequencies, sd = stds, scenario = scenarios, jumpSize = jumpSizes)




##### FUNCTION FOR RUNNING SIMULATIONS #####

runSim <- function(i, simulations) {
  fileName <- paste(c("simulations/resSine/", simulations[i, ], ".RData"), collapse = "")
  if (!file.exists(fileName)) {
    cat("Running ", fileName, "\n")
    p <- simulations[i, ]

    Y <- mclapply(1:REPS, function(r) dataSinusoidal(N, amplitude = p$amplitude, frequency = p$frequency, sd = p$sd, type = as.character(p$scenario), jumpSize = p$jumpSize), mc.cores = CORES)
    signal <- lapply(Y, function(r) r$signal)
    y <- lapply(Y, function(r) r$y)
    changepoints <- Y[[1]]$changepoints

    # DeCAFS K 15
    resDeCAFSESTK15 <- mclapply(y, DeCAFS, mc.cores = CORES)
    
    # ar1seg with estimator
    resar1segEST <- mclapply(y, AR1seg_func, Kmax = 10, mc.cores = CORES)
    
    # threshold with estimator
    resThresholdEST15 <- lapply(y, function(y){
      est <- estimateParameters(y)
      l2Threshold(y, beta = 2 * log(N),  lambda =  1/est$sdEta^2)
    })
    
    
    save(signal,
         y,
         changepoints,
         resDeCAFSESTK15,
         resar1segEST,
         resThresholdEST15,
         file = fileName)
  }
}



##### FIXED AMPLITUDES, RANGING FREQUENCIES #####

toSummarize <- simulations %>% filter(amplitude == 15, sd == 2)
if (T) lapply(1:nrow(toSummarize), runSim, simulations = toSummarize)


# summary df
F1df <- lapply(1:nrow(toSummarize), function(i) {
  p <- toSummarize[i, ]
  print(p)
  
  fileName <- paste(c("simulations/resSine/", p, ".RData"), collapse = "")
  if (!file.exists(fileName)) {
    cat("Missing", paste0(Map(paste, names(p), p), collapse = " "), "\n")
    return(NULL)
  } else load(fileName)
  
  DeCAFSdfK15 <- cbind(p$frequency,
                       sapply(resDeCAFSESTK15, function(r)
                         computeF1Score(c(changepoints, N), c(r$changepoints,N), 3)) %>% as.numeric,
                       sapply(resDeCAFSESTK15, function(r)
                         computePrecision(c(changepoints, N), c(r$changepoints,N), 3)) %>% as.numeric,
                       sapply(resDeCAFSESTK15, function(r)
                         computeRecall(c(changepoints, N), c(r$changepoints,N), 3)) %>% as.numeric,
                       as.character(p$scenario), 
                       "DeCAFS est")
  AR1segdf <- cbind(p$frequency,
                     sapply(resar1segEST, function(r)
                       computeF1Score(c(changepoints, N), r$PPSelectedBreaks, 3)) %>% as.numeric,
                     sapply(resar1segEST, function(r)
                       computePrecision(c(changepoints, N), r$PPSelectedBreaks, 3)) %>% as.numeric, # comptuting the F1
                     sapply(resar1segEST, function(r)
                       computeRecall(c(changepoints, N), r$PPSelectedBreaks, 3)) %>% as.numeric, # comptuting the F1
                    as.character(p$scenario), 
                    "AR1Seg est")
  return(rbind(AR1segdf, DeCAFSdfK15))
})


F1df <- Reduce(rbind, F1df)

colnames(F1df) <- c("frequency", "F1Score", "Precision", "Recall", "Scenario", "Algorithm")
F1df <- as_tibble(F1df) %>% mutate(frequency = as.numeric(frequency),
                                   F1Score = as.numeric(F1Score),
                                   Precision = as.numeric(Precision),
                                    Recall = as.numeric(Recall))


save(F1df, file = "simulations/outputs/F1Sinusoidal.RData")
load("simulations/outputs/F1Sinusoidal.RData")


cbPalette3 <- c("#56B4E9",  "#33cc00")
scores <- ggplot(F1df, aes(x = frequency, y = F1Score, group = Algorithm, color = Algorithm, by = Algorithm)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .0001) +
  facet_wrap(. ~ Scenario ) + 
  scale_color_manual(values = cbPalette3)
scores




# example plot
p <- toSummarize[12, ]
fileName <- paste(c("simulations/resSine/", p, ".RData"), collapse = "")
if (!file.exists(fileName)) {
  cat("Missing", paste0(Map(paste, names(p), p), collapse = " "), "\n")
  return(NULL)
} else load(fileName)



k <- 48

# estimated spikes DeCAFS
df <- data.frame(x1 = resDeCAFSESTK15[[k]]$changepoints, y1 = -20, y2 = -25)
estimDeCAFS = geom_segment(aes(x = x1, xend = x1, y = y1, yend = y2), data = df, col = cbPalette3[2])


# estimated spikes AR(1)Seg
df <- data.frame(x1 = resar1segEST[[k]]$PPSelectedBreaks[1:(length(resar1segEST[[k]]$PPSelectedBreaks) - 1)], y1 = -25, y2 = -30)
estimAR1Seg = geom_segment(aes(x = x1, xend = x1, y = y1, yend = y2), data = df, col = cbPalette3[1])



exe <- ggplot(data.frame(t = 1:length(y[[k]]), y[[k]]), aes(x = t, y = y[[k]])) +
  geom_point(col = "grey") +
  geom_line(aes(x = t, y = value, lty = signal), data = data.frame(t = 1:length(y[[k]]), real = signal[[k]]) %>% gather(signal, value, -t)) +
  ylab("y") + 
  theme(legend.position = "none")

exe + estimDeCAFS + estimAR1Seg

compositePlot <- ggarrange(scores, exe + estimDeCAFS + estimAR1Seg, heights = c(2,1), labels = c("A", "B"), ncol = 1, legend = "top",
    common.legend = T)

ggsave(compositePlot, width = 9, height = 10, file = "simulations/outputs/4-Sinusoidal.pdf", device = "pdf", dpi = "print")



# adding plots for precision and recall
Prec <- ggplot(F1df, aes(x = frequency, y = Precision, group = Algorithm, color = Algorithm, by = Algorithm)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .0001) +
  facet_wrap(. ~ Scenario ) +
  scale_color_manual(values = cbPalette3)
Prec
ggsave(Prec, width = 6, height = 4, units = "in", file = "simulations/outputs/4-SinusoidalPrec.pdf", device = "pdf", dpi = "print")

# adding plots for precision and recall
Rec <- ggplot(F1df, aes(x = frequency, y = Recall, group = Algorithm, color = Algorithm, by = Algorithm)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .0001) +
  facet_wrap(. ~ Scenario ) +
  scale_color_manual(values = cbPalette3)
Rec
ggsave(Rec, width = 6, height = 4, units = "in", file = "simulations/outputs/4-SinusoidalRec.pdf", device = "pdf", dpi = "print")
