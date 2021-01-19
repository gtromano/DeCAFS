# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                           #
#     Return the plot for the precision and recall from the previous        #
#     simulations - REQUIRES 1-RWAR FIRST                                   #
#                                                                           #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
library(ggpubr)

# ~~~~~~~~~~~~~~~~ #
#    precision     #
# ~~~~~~~~~~~~~~~~ #


# first batch
toSummarize <- simulations %>% filter(sigmaEta == 0, sigmaNu == 2, jumpSize == 10)
F1prec <- mclapply(1:nrow(toSummarize), function(i) {
  p <- toSummarize[i, ]
  print(p)

  fileName <- paste(c("simulations/resRWAR/", p, ".RData"), collapse = "")
  if (!file.exists(fileName)) {
    cat("Missing", paste0(Map(paste, names(p), p), collapse = " "), "\n")
    return(NULL)
  } else load(fileName)

  DeCAFSdf <- cbind(p$phi,
                    sapply(resDeCAFS, function(r) computePrecision(c(changepoints, N), c(r$changepoints, N), 3)) %>% as.numeric,
                    as.character(p$scenario),
                    "DeCAFS")
  fpopdf <- cbind(p$phi,
                  sapply(resfpop, function(r) computePrecision(c(changepoints, N), r$t.est, 3)) %>% as.numeric,
                  as.character(p$scenario),
                  "fpop")
  enffpopdf <- cbind(p$phi,
                     sapply(resenffpop, function(r) computePrecision(c(changepoints, N), r$t.est, 3)) %>% as.numeric,
                     as.character(p$scenario),
                     "fpop Inf")
  AR1segdf <- cbind(p$phi,
                    sapply(resar1seg, function(r)
                      computePrecision(c(changepoints, N), r$PPSelectedBreaks, 3)) %>% as.numeric, # comptuting the F1
                    as.character(p$scenario),
                    "AR1Seg")

  AR1segdfest <- cbind(p$phi,
                       sapply(resar1segEST, function(r)
                         computePrecision(c(changepoints, N), r$PPSelectedBreaks, 3)) %>% as.numeric, # comptuting the F1
                       as.character(p$scenario),
                       "AR1Seg est")
  DeCAFSdfK15 <- cbind(p$phi,
                       sapply(resDeCAFSESTK15, function(r) computePrecision(c(changepoints, N), c(r$changepoints, N), 3)) %>% as.numeric,
                       as.character(p$scenario),
                       "DeCAFS est")

  NPPELT <- cbind(p$phi,
                  sapply(resNPPELT, function(r) computePrecision(c(changepoints, N), c(r@cpts, N), 3)) %>% as.numeric,
                  as.character(p$scenario),
                  "NP-PELT")

  return(rbind(DeCAFSdf, fpopdf, enffpopdf, AR1segdf, DeCAFSdfK15, AR1segdfest, NPPELT))
}, mc.cores = 6)

F1prec <- Reduce(rbind, F1prec)
colnames(F1prec) <- c("phi", "Precision", "Scenario", "Algorithm")
F1prec <- as_tibble(F1prec) %>% mutate(phi = as.numeric(phi),
                                   Precision = as.numeric(Precision))

cbPalette <- c("#0072B2", "#56B4E9", "#009E73", "#33cc00", "#E69F00", "#CC79A7", "#984447")
scores <- ggplot(F1prec %>% filter(Algorithm != "DeCAFS est (5)" & Algorithm != "DeCAFS est (10)"),
                 aes(x = phi, y = Precision, group = Algorithm, by = Algorithm, col = Algorithm)) +
  geom_vline(xintercept = unique(simulations$phi)[7], col = "grey", lty = 2) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .05) +
  facet_wrap(~ Scenario) +
  scale_color_manual(values = cbPalette) +
  xlab(expression(italic(phi)))
scores



# second batch
toSummarize <- simulations %>% filter(sigmaEta == 0, sigmaNu == 2, jumpSize == 10)
F1prec <- lapply(1:nrow(toSummarize), function(i) {
  p <- toSummarize[i, ]
  print(p)

  fileName <- paste(c("simulations/resRWAR/", p, ".RData"), collapse = "")
  if (!file.exists(fileName)) {
    cat("Missing", paste0(Map(paste, names(p), p), collapse = " "), "\n")
    return(NULL)
  } else load(fileName)

  DeCAFSdf <- cbind(p$jumpSize,
                    sapply(resDeCAFS, function(r) computePrecision(c(changepoints, N), c(r$changepoints, N), 3)) %>% as.numeric,
                    as.character(p$scenario),
                    "DeCAFS")
  fpopdf <- cbind(p$jumpSize,
                  sapply(resfpop, function(r) computePrecision(c(changepoints, N), r$t.est, 3)) %>% as.numeric,
                  as.character(p$scenario),
                  "fpop")
  enffpopdf <- cbind(p$jumpSize,
                     sapply(resenffpop, function(r) computePrecision(c(changepoints, N), r$t.est, 3)) %>% as.numeric,
                     as.character(p$scenario),
                     "fpop Inf")
  AR1segdf <- cbind(p$jumpSize,
                    sapply(resar1seg, function(r)
                      computePrecision(c(changepoints, N), r$PPSelectedBreaks, 3)) %>% as.numeric, # comptuting the F1
                    as.character(p$scenario),
                    "AR1Seg")
  AR1segdfest <- cbind(p$jumpSize,
                       sapply(resar1segEST, function(r)
                         computePrecision(c(changepoints, N), r$PPSelectedBreaks, 3)) %>% as.numeric, # comptuting the F1
                       as.character(p$scenario),
                       "AR1Seg est")
  DeCAFSdfK15 <- cbind(p$jumpSize,
                       sapply(resDeCAFSESTK15, function(r) computePrecision(c(changepoints, N), c(r$changepoints, N), 3)) %>% as.numeric,
                       as.character(p$scenario),
                       "DeCAFS est")

  NPPELT <- cbind(p$jumpSize,
                sapply(resNPPELT, function(r) computePrecision(c(changepoints, N), c(r@cpts, N), 3)) %>% as.numeric,
                as.character(p$scenario),
                "NP-PELT")

  return(rbind(DeCAFSdf, fpopdf, enffpopdf, AR1segdf,  DeCAFSdfK15, AR1segdfest, NPPELT))
})


F1prec <- Reduce(rbind, F1prec)

colnames(F1prec) <- c("jumpSize", "Precision", "Scenario", "Algorithm")
F1prec <- as_tibble(F1prec) %>% mutate(jumpSize = as.numeric(jumpSize),
                                   Precision = as.numeric(Precision))
cbPalette <- c("#0072B2", "#56B4E9", "#009E73", "#33cc00", "#E69F00", "#CC79A7", "#984447")
scoresJump <- ggplot(F1prec %>% filter(Algorithm != "DeCAFS est (5)" & Algorithm != "DeCAFS est (10)"),
                     aes(x = jumpSize, y = Precision, group = Algorithm, by = Algorithm, col = Algorithm)) +
  geom_vline(xintercept = 10, col = "grey", lty = 2) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .5) +
  facet_wrap(~ Scenario) +
  scale_color_manual(values = cbPalette) +
  xlab("Jump Size")
scoresJump



# third batch

# selecting the simulations with sd_nu = 2
toSummarize <- simulations %>% filter(phi == (simulations$phi %>% unique())[7], sigmaNu == 2, jumpSize == 10)

if (F) mclapply(1:nrow(toSummarize), runSim, simulations = toSummarize, SIMID = "", mc.cores = 8)


F1prec <- mclapply(1:nrow(toSummarize), function(i) {
  p <- toSummarize[i, ]
  print(p)

  fileName <- paste(c("simulations/resRWAR/", p, ".RData"), collapse = "")
  if (!file.exists(fileName)) {
    cat("Missing", paste0(Map(paste, names(p), p), collapse = " "), "\n")
    return(NULL)
  } else load(fileName)

  DeCAFSdf <- cbind(p$sigmaEta,
                    sapply(resDeCAFS, function(r) computePrecision(c(changepoints, N), c(r$changepoints, N), 3)) %>% as.numeric,
                    as.character(p$scenario),
                    "DeCAFS")
  fpopdf <- cbind(p$sigmaEta,
                  sapply(resfpop, function(r) computePrecision(c(changepoints, N), r$t.est, 3)) %>% as.numeric,
                  as.character(p$scenario),
                  "fpop")
  enffpopdf <- cbind(p$sigmaEta,
                     sapply(resenffpop, function(r) computePrecision(c(changepoints, N), r$t.est, 3)) %>% as.numeric,
                     as.character(p$scenario),
                     "fpop Enf")
  AR1segdf <- cbind(p$sigmaEta,
                    sapply(resar1seg, function(r)
                      computePrecision(c(changepoints, N), r$PPSelectedBreaks, 3)) %>% as.numeric, # comptuting the F1
                    as.character(p$scenario),
                    "AR1Seg")
  AR1segdfest <- cbind(p$sigmaEta,
                       sapply(resar1segEST, function(r)
                         computePrecision(c(changepoints, N), r$PPSelectedBreaks, 3)) %>% as.numeric, # comptuting the F1
                       as.character(p$scenario),
                       "AR1Seg est")
  DeCAFSdfK15 <- cbind(p$sigmaEta,
                       sapply(resDeCAFSESTK15, function(r) computePrecision(c(changepoints, N), c(r$changepoints, N), 3)) %>% as.numeric,
                       as.character(p$scenario),
                       "DeCAFS est")

  NPPELT <- cbind(p$sigmaEta,
                sapply(resNPPELT, function(r) computePrecision(c(changepoints, N), c(r@cpts, N), 3)) %>% as.numeric,
                as.character(p$scenario),
                "NP-PELT")

  return(rbind(DeCAFSdf, fpopdf, enffpopdf, AR1segdf, DeCAFSdfK15, AR1segdfest, NPPELT))
}, mc.cores = 4)


F1prec <- Reduce(rbind, F1prec)

colnames(F1prec) <- c("SigmaEta", "Precision", "Scenario", "Algorithm")
F1prec <- as_tibble(F1prec) %>% mutate(SigmaEta = as.numeric(SigmaEta),
                                   Precision = as.numeric(Precision))

cbPaletteEDIT <- c("#56B4E9", "#009E73", "#33cc00", "#E69F00", "#CC79A7", "#984447")
scoresRWAR <- ggplot(F1prec %>% filter(Algorithm != "DeCAFS est (5)" & Algorithm != "DeCAFS est (10)" & Algorithm != "AR1Seg"),
                     aes(x = SigmaEta, y = Precision, group = Algorithm, by = Algorithm, col = Algorithm)) +
  geom_vline(xintercept = 0, col = "grey", lty = 2) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.08) +
  facet_wrap(~ Scenario) +
  scale_color_manual(values = cbPaletteEDIT) +
  xlab(expression(sigma[eta]))

scoresRWAR
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
