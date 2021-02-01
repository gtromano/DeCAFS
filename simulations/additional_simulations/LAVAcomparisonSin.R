# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                       #
#   The LAVA comparison on a sinusoidal process with abrupt changes     #
#                                                                       #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

library(tidyverse)
library(parallel) # for mclapply
library(DeCAFS)
library(Lavash)
library(AR1seg)   # for comparison with AR1seg

source("simulations/helper_functions.R")

lavaCHANGEPOINT <- function(y, l1penalty, l2penalty) {
  N <- length(y)
  # creating a lower triangular matrix
  L <- matrix(1, N, N)
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


REPS <- 100 # number of replicates
N <- 1e3 # lenght of the sequence
CORES <- 16

# range of model parameters
amplitudes <- seq(1, 5, length.out = 5)
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
  fileName <- paste(c("simulations/additional_simulations/resLAVA/", simulations[i, ], ".RData"), collapse = "")
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

    # DeCAFS with tweaked penalties
    resDeCAFS <- mclapply(y, DeCAFS, beta = (2 * log(N)), modelParam = list(sdEta = estVariation, sdNu = p$sd, phi = 0), mc.cores = CORES)

    # DeCAFS K 15
    resDeCAFSESTK15 <- lapply(y, DeCAFS)

    # LAVA oracle
    resLAVA <- mclapply(y, lavaCHANGEPOINT, l1penalty = seq(.1, 1, length.out = 40), l2penalty = getLavaPenalty(estVariation, p$sd, N), mc.cores = CORES)

    # LAVA with the same estimates as DeCAFS
    params <- lapply (y, estimateParameters)
    resLAVAESTK15 <- mclapply(1:REPS, function (r) lavaCHANGEPOINT(y[[r]], l1penalty = seq(.1, 1, length.out = 40), l2penalty = getLavaPenalty(params[[r]]$sdEta, params[[r]]$sdNu, N)), mc.cores = CORES)

    save(signal,
         y,
         changepoints,
         resDeCAFS,
         resDeCAFSESTK15,
         resLAVA,
         resLAVAESTK15,
         file = fileName)
  }
}


# running simulations
toSummarize <- simulations %>% filter(nbSeg == 20)


if (T) lapply(1:nrow(toSummarize), runSim, simulations = toSummarize)


# summary df
df <- lapply(1:nrow(toSummarize), function(i) {
  p <- toSummarize[i, ]
  print(p)

  fileName <- paste(c("simulations/additional_simulations/resLAVA/", p, ".RData"), collapse = "")
  if (!file.exists(fileName)) {
    cat("Missing", paste0(Map(paste, names(p), p), collapse = " "), "\n")
    return(NULL)
  } else load(fileName)

  DeCAFSdf <- cbind(p$amplitude,
                     sapply(resDeCAFS, function(r)
                       computeF1Score(c(changepoints, N), c(r$changepoints,N), 3)) %>% as.numeric,
                     sapply(resDeCAFS, function(r)
                       computePrecision(c(changepoints, N), c(r$changepoints,N), 3)) %>% as.numeric,
                     sapply(resDeCAFS, function(r)
                       computeRecall(c(changepoints, N), c(r$changepoints,N), 3)) %>% as.numeric,
                     sapply(resDeCAFS, function (r) {
                       mse(signal[[1]], r$signal)
                      }),
                     as.character(p$scenario),
                     "DeCAFS")


  DeCAFSdfK15 <- cbind(p$amplitude,
                       sapply(resDeCAFSESTK15, function(r)
                         computeF1Score(c(changepoints, N), c(r$changepoints,N), 3)) %>% as.numeric,
                       sapply(resDeCAFSESTK15, function(r)
                         computePrecision(c(changepoints, N), c(r$changepoints,N), 3)) %>% as.numeric,
                       sapply(resDeCAFSESTK15, function(r)
                         computeRecall(c(changepoints, N), c(r$changepoints,N), 3)) %>% as.numeric,
                       sapply(resDeCAFSESTK15, function (r) {
                         mse(signal[[1]], r$signal)
                        }),
                       as.character(p$scenario),
                       "DeCAFS est")

  LAVAdf <- cbind(p$amplitude,
                    sapply(resLAVA, function(r)
                       computeF1Score(c(changepoints, N), c(r$cp, N), 3)) %>% as.numeric,
                    sapply(resLAVA, function(r)
                       computePrecision(c(changepoints, N), c(r$cp, N), 3)) %>% as.numeric,
                    sapply(resLAVA, function(r)
                       computeRecall(c(changepoints, N), c(r$cp, N), 3)) %>% as.numeric,
                    sapply(resLAVA, function (r) {
                         mse(signal[[1]], r$est)
                       }),
                    as.character(p$scenario),
                      "LAVA")

  LAVAdfest <- cbind(p$amplitude,
                  sapply(resLAVAESTK15, function(r)
                     computeF1Score(c(changepoints, N), c(r$cp, N), 3)) %>% as.numeric,
                  sapply(resLAVAESTK15, function(r)
                     computePrecision(c(changepoints, N), c(r$cp, N), 3)) %>% as.numeric,
                  sapply(resLAVAESTK15, function(r)
                     computeRecall(c(changepoints, N), c(r$cp, N), 3)) %>% as.numeric,
                  sapply(resLAVAESTK15, function (r) {
                       mse(signal[[1]], r$est)
                     }),
                  as.character(p$scenario),
                    "LAVA est")


  return(rbind(DeCAFSdf, DeCAFSdfK15, LAVAdf, LAVAdfest))
})


df <- Reduce(rbind, df)

colnames(df) <- c("amplitude", "F1Score", "Precision", "Recall", "mse", "Scenario", "Algorithm")
df <- as_tibble(df) %>% mutate(amplitude = as.numeric(amplitude),
                               F1Score = as.numeric(F1Score),
                               Precision = as.numeric(Precision),
                               Recall = as.numeric(Recall),
                               mse = as.numeric(mse))


#save(df, file = "simulations/additional_simulations/resLAVA/df.RData")
#load("simulations/additional_simulations/resLAVA/df.RData")


cbPalette3 <- c("#009E73", "#33cc00", "#ef0716", "#fc7658")
cbPalette4 <- c("#33cc00", "#ef0716", "#fc7658")

F1df <- ggplot(df %>% filter(Algorithm != "DeCAFS"), aes(x = amplitude, y = F1Score, group = Algorithm, color = Algorithm, by = Algorithm)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .001) +
  facet_wrap(. ~ Scenario ) +
  ylim(0, 1) +
  scale_color_manual(values = cbPalette4)

Prec <- ggplot(df %>% filter(Algorithm != "DeCAFS"), aes(x = amplitude, y = Precision, group = Algorithm, color = Algorithm, by = Algorithm)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .001) +
  facet_wrap(. ~ Scenario ) +
  ylim(0, 1) +
  scale_color_manual(values = cbPalette4)

Recall <- ggplot(df %>% filter(Algorithm != "DeCAFS"), aes(x = amplitude, y = Recall, group = Algorithm, color = Algorithm, by = Algorithm)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .001) +
  facet_wrap(. ~ Scenario ) +
  ylim(0, 1) +
  scale_color_manual(values = cbPalette4)


mse <- ggplot(df %>% filter(Algorithm != "DeCAFS"), aes(x = amplitude, y = mse, group = Algorithm, color = Algorithm, by = Algorithm)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .001) +
  ylab("MSE") +
  facet_wrap(. ~ Scenario ) +
  scale_color_manual(values = cbPalette4)



# # example plot
# p <- toSummarize[2, ]
# fileName <- paste(c("simulations/additional_simulations/resLAVA/", p, ".RData"), collapse = "")
# if (!file.exists(fileName)) {
#   cat("Missing", paste0(Map(paste, names(p), p), collapse = " "), "\n")
#   return(NULL)
# } else load(fileName)
#
# k <- 42
# # estimated spikes DeCAFS
# df2 <- data.frame(x1 = resDeCAFSESTK15[[k]]$changepoints, y1 = -10, y2 = -15)
# estimDeCAFS = geom_segment(aes(x = x1, xend = x1, y = y1, yend = y2), data = df2, col = cbPalette4[1])
#
# # estimated spikes AR(1)Seg
# df2 <- data.frame(x1 = resLAVA[[k]]$cp, y1 = -15, y2 = -20)
# estimAR1Seg = geom_segment(aes(x = x1, xend = x1, y = y1, yend = y2), data = df2, col = cbPalette4[2])
#
# y1 <- y[[k]]
# exe <- ggplot(data.frame(t = 1:length(y[[k]]), y1), aes(x = t, y = y1)) +
#   geom_point(col = "grey") +
#   geom_line(aes(x = t, y = value, color = signal), data = data.frame(t = 1:length(y[[k]]), DeCAFS = resDeCAFSESTK15[[k]]$signal, LAVA = resLAVA[[k]]$est, LAVAEST = resLAVAESTK15[[k]]$est) %>% gather(signal, value, -t)) +
#   ylab("y") +
#   scale_color_manual(values = cbPalette4) +
#   xlim(0, 250) +
#   ylim(-20, 15) +
#   theme(legend.position = "none")
#
# example1 <- exe + estimDeCAFS + estimAR1Seg
#
#
### example plot 2
p <- toSummarize[2, ]
fileName <- paste(c("simulations/additional_simulations/resLAVA/", p, ".RData"), collapse = "")
if (!file.exists(fileName)) {
  cat("Missing", paste0(Map(paste, names(p), p), collapse = " "), "\n")
  return(NULL)
} else load(fileName)

k <- 42
# estimated spikes DeCAFS
df2 <- data.frame(x1 = resDeCAFSESTK15[[k]]$changepoints, y1 = -10, y2 = -15)
estimDeCAFS = geom_segment(aes(x = x1, xend = x1, y = y1, yend = y2), data = df2, col = cbPalette4[1])

# estimated spikes AR(1)Seg
df2 <- data.frame(x1 = resLAVA[[k]]$cp, y1 = -15, y2 = -20)
estimAR1Seg = geom_segment(aes(x = x1, xend = x1, y = y1, yend = y2), data = df2, col = cbPalette4[2])

y2 <- y[[k]]
exe <- ggplot(data.frame(t = 1:length(y[[k]]), y[[k]]), aes(x = t, y = y2)) +
  geom_point(col = "grey") +
  geom_line(aes(x = t, y = value, color = signal, alpha = .8), data = data.frame(t = 1:length(y[[k]]), DeCAFS = resDeCAFSESTK15[[k]]$signal, LAVA = resLAVA[[k]]$est, LAVAEST = resLAVAESTK15[[k]]$est) %>% gather(signal, value, -t)) +
  ylab("y") +
  scale_color_manual(values = cbPalette4) +
  xlim(0, 250) +
  ylim(-20, 15) +
  theme(legend.position = "none")

example2 <- exe + estimDeCAFS + estimAR1Seg


### composite plot
library(ggpubr)

# getting legend out

meaplot <- ggarrange(
    F1df,
    Prec,
    Recall,
    labels = c("A1", "A2", "A3"),
    ncol = 3,
    legend = "top",
    common.legend = T
  )


exeplot <- ggarrange(
  example1 + xlim(0, 250), example2,
  labels = c("B1", "B2")
)


ggsave(ggarrange(meaplot, exeplot, ncol = 1), width = 9, height = 8, file = "simulations/outputs/LAVAcomp.pdf", device = "pdf", dpi = "print")
