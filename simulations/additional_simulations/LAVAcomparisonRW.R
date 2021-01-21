# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                               #
#   The LAVA comparison on a RW process with abrupt changes     #
#                                                               #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

library(tidyverse)
library(parallel) # for mclapply
library(DeCAFS)
library(Lavash)

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


getLavaPenalty <- function (sdEta, sdNu, N) sdNu^2 * (1 / N * sdEta ^ 2)


REPS <- 100 # number of replicates
N <- 1e3 # lenght of the sequence
CORES <- 16

# generate a list of simulations
simulations <- expand.grid(sdEta = seq(0, 2, length.out = 5), sdNu = 2,  jumpSize = 15, scenario = "updown", nbSeg = 20)


##### FUNCTION FOR RUNNING SIMULATIONS #####

runSim <- function(i, simulations) {
  fileName <- paste(c("simulations/additional_simulations/resLAVARW/", simulations[i, ], ".RData"), collapse = "")
  if (!file.exists(fileName)) {

    cat("Running ", fileName, "\n")
    p <- simulations[i, ]
    Y <- mclapply(1:REPS, function(r) dataRWAR(N, sdEta = p$sdEta, sdNu = p$sdNu,  jumpSize = p$jumpSize, type = as.character(p$scenario), nbSeg = p$nbSeg), mc.cores = CORES)

    signal <- lapply(Y, function(r) r$signal)
    y <- lapply(Y, function(r) r$y)
    changepoints <- Y[[1]]$changepoints

    resDeCAFS <- mclapply(y, DeCAFS, beta = (2 * log(N)), modelParam = list(sdEta = p$sdEta, sdNu = p$sdNu, phi = 0), mc.cores = CORES)

    # DeCAFS K 15 on Random Walk
    params <- lapply (y, function(y_s) {
        pest <- estimateParameters(y_s, phiUpper = 1e-10)
        pest$phi <- 0
        return(pest)
      }
    )
    resDeCAFSESTK15 <- mclapply(1:REPS, function(r) DeCAFS(y[[r]], modelParam = params[[r]]), mc.cores = CORES)

    # LAVA oracle
    resLAVA <- mclapply(y, lavaCHANGEPOINT, l1penalty = seq(.1, 1, length.out = 20), l2penalty = getLavaPenalty(p$sdEta, p$sdNu, N), mc.cores = CORES)

    # LAVA with the same estimates as DeCAFS
    resLAVAESTK15 <- mclapply(1:REPS, function (r) lavaCHANGEPOINT(y[[r]], l1penalty = seq(.1, 1, length.out = 20), l2penalty = getLavaPenalty(params[[r]]$sdEta, params[[r]]$sdNu, N)), mc.cores = CORES)

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


### running simulations ####
toSummarize <- simulations %>% filter(nbSeg == 20)


if (T) lapply(1:nrow(toSummarize), runSim, simulations = toSummarize)


#### SUMMARIZING SIMULATIONS ####
# summary df
df <- lapply(1:nrow(toSummarize), function(i) {
  p <- toSummarize[i, ]
  print(p)

  fileName <- paste(c("simulations/additional_simulations/resLAVARW/", p, ".RData"), collapse = "")
  if (!file.exists(fileName)) {
    cat("Missing", paste0(Map(paste, names(p), p), collapse = " "), "\n")
    return(NULL)
  } else load(fileName)

  DeCAFS <- cbind(p$sdEta,
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


  DeCAFSdfK15 <- cbind(p$sdEta,
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

  LAVAdf <- cbind(p$sdEta,
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

    LAVAdfest <- cbind(p$sdEta,
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

  return(rbind(DeCAFS, DeCAFSdfK15, LAVAdf, LAVAdfest))
})


df <- Reduce(rbind, df)

colnames(df) <- c("sigmaEta", "F1Score", "Precision", "Recall", "mse", "Scenario", "Algorithm")
df <- as_tibble(df) %>% mutate(sigmaEta = as.numeric(sigmaEta),
                               F1Score = as.numeric(F1Score),
                               Precision = as.numeric(Precision),
                               Recall = as.numeric(Recall),
                               mse = as.numeric(mse))


# save(df, file = "simulations/additional_simulations/resLAVARW/dfRW.RData")
# load("simulations/additional_simulations/resLAVARW/dfRW.RData")


cbPalette3 <- c("#009E73", "#33cc00", "#56B4E9", "#0072B2")
Prec <- ggplot(df, aes(x = sigmaEta, y = Precision, group = Algorithm, color = Algorithm, by = Algorithm)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .001) +
  facet_wrap(. ~ Scenario ) +
  ylim(0, 1) +
  scale_color_manual(values = cbPalette3) +
  xlab(expression(sigma[eta]))


Recall <- ggplot(df, aes(x = sigmaEta, y = Recall, group = Algorithm, color = Algorithm, by = Algorithm)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .001) +
  facet_wrap(. ~ Scenario ) +
  ylim(0, 1) +
  scale_color_manual(values = cbPalette3) +
  xlab(expression(sigma[eta]))


mse <- ggplot(df, aes(x = sigmaEta, y = mse, group = Algorithm, color = Algorithm, by = Algorithm)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .001) +
  ylab("MSE") +
  facet_wrap(. ~ Scenario ) +
  scale_color_manual(values = cbPalette3) +
  xlab(expression(sigma[eta]))



# example plot
p <- toSummarize[2, ]
fileName <- paste(c("simulations/additional_simulations/resLAVARW/", p, ".RData"), collapse = "")
if (!file.exists(fileName)) {
  cat("Missing", paste0(Map(paste, names(p), p), collapse = " "), "\n")
  return(NULL)
} else load(fileName)

k <- 42
# estimated spikes DeCAFS
df2 <- data.frame(x1 = resDeCAFS[[k]]$changepoints, y1 = -10, y2 = -15)
estimDeCAFS = geom_segment(aes(x = x1, xend = x1, y = y1, yend = y2), data = df2, col = cbPalette3[1])

# estimated spikes AR(1)Seg
df2 <- data.frame(x1 = resLAVA[[k]]$cp, y1 = -15, y2 = -20)
estimAR1Seg = geom_segment(aes(x = x1, xend = x1, y = y1, yend = y2), data = df2, col = cbPalette3[3])

y1 <- y[[k]]
exe <- ggplot(data.frame(t = 1:length(y[[k]]), y[[k]]), aes(x = t, y = y1)) +
  geom_point(col = "grey") +
  geom_line(aes(x = t, y = value, color = signal), data = data.frame(t = 1:length(y[[k]]), DeCAFS = resDeCAFSESTK15[[k]]$signal, LAVA = resLAVA[[k]]$est) %>% gather(signal, value, -t)) +
  ylab("y") +
  scale_color_manual(values = cbPalette3) +
  xlim(0, 250) +
  theme(legend.position = "none")

example1 <- exe + estimDeCAFS + estimAR1Seg


### example plot 2
p <- toSummarize[5, ]
fileName <- paste(c("simulations/additional_simulations/resLAVARW/", p, ".RData"), collapse = "")
if (!file.exists(fileName)) {
  cat("Missing", paste0(Map(paste, names(p), p), collapse = " "), "\n")
  return(NULL)
} else load(fileName)

k <- 42
# estimated spikes DeCAFS
df2 <- data.frame(x1 = resDeCAFSESTK15[[k]]$changepoints, y1 = -20, y2 = -40)
estimDeCAFS = geom_segment(aes(x = x1, xend = x1, y = y1, yend = y2), data = df2, col = cbPalette3[1])

# estimated spikes AR(1)Seg
df2 <- data.frame(x1 = resLAVA[[k]]$cp, y1 = -40, y2 = -60)
estimAR1Seg = geom_segment(aes(x = x1, xend = x1, y = y1, yend = y2), data = df2, col = cbPalette3[3])

y2 <- y[[k]]
exe <- ggplot(data.frame(t = 1:length(y[[k]]), y[[k]]), aes(x = t, y = y2)) +
  geom_point(col = "grey") +
  geom_line(aes(x = t, y = value, color = signal), data = data.frame(t = 1:length(y[[k]]), DeCAFS = resDeCAFSESTK15[[k]]$signal, LAVA = resLAVA[[k]]$est) %>% gather(signal, value, -t)) +
  ylab("y") +
  scale_color_manual(values = cbPalette3) +
  xlim(0, 250) +
  ylim(-60, 65) +
  theme(legend.position = "none")

example2 <- exe + estimDeCAFS + estimAR1Seg

### composite plot
library(ggpubr)

# getting legend out

meaplot <- ggarrange(
    Prec,
    Recall,
    mse + scale_y_log10() + ylab("mse (log scale)"),
    labels = c("A1", "A2", "A3"),
    ncol = 3,
    legend = "top",
    common.legend = T
  )


exeplot <- ggarrange(
  example1, example2,
  labels = c("B1", "B2")
)


ggarrange(meaplot, exeplot, ncol = 1)

ggsave(ggarrange(meaplot, exeplot, ncol = 1), width = 8, height = 7, file = "simulations/outputs/LAVAcompRW.pdf", device = "pdf", dpi = "print")

