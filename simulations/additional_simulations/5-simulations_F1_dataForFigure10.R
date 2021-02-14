rm(list=ls())

# LIBRARIES
#devtools::install_github("gtromano/DeCAFS")
#devtools::install_github("vrunge/ARRWestim")
library(parallel)
library(fields)

###
### NB cores for parallel computing
###
cores <- detectCores()
cores <- 60
#cores <- 8 ###CHANGE


one.simu <- function(i, N = 10^5, sdEta = 0.4, sdNu = 0.3, phi = 0.2,
                     type = "rand1", nbSeg = 10, jumpSize = 2, nbK = 10, varType = "MAD")
{

  y <- ARRWestim::dataRWAR(N = N,
                sdEta = sdEta, sdNu = sdNu, phi = phi,
                type = type,
                nbSeg = nbSeg, jumpSize = jumpSize,
                seed = sample(1e6,1))
  res <- bestParameters(y$y, nbK = nbK, type = varType)

  #ind, n, phi, sdEta2, sdNu2, nbK, poiP, meanGap, method, mse, beta, K
  df <- data.frame(numeric(0), numeric(0), numeric(0), numeric(0), numeric(0),
                   numeric(0), numeric(0), numeric(0), numeric(0), numeric(0),
                   numeric(0) ,stringsAsFactors = FALSE)
  colnames(df) <- c("index", "n", "sdEta", "sdNu", "phi", "nbK",
                    "nbSeg", "jumpSize", "sdEtaEst%", "sdNuEst%", "phiEst_error")
  df[1,] <- c(i,N, sdEta, sdNu, phi, nbK, nbSeg, jumpSize,
              (res$EtaOpt - sdEta)/sdEta, (res$NuOpt-sdNu)/sdNu, res$argmin - phi)
  return(df)
}

########### ########### ########### ###########

###
### Simulations parameters
###
#nbSimu <- 1200
nbSimu <- 2000 ##CHANGE
nbPhi <- 18 #step size = 0.05 in phi
nbOmega2 <- 40
nbK <- 10


###
### phi and omega2 = grid for simulations (omega2 = (sd_eta/sd_nu)^2)
### sd_nu fixed to 1
###

phi <- seq(from = 0, to = 0.85, length.out = nbPhi)

omega2 <- exp(seq(from = -log(12), to = log(2), length.out = nbOmega2))
diffO2 <- diff(log(omega2))[1]
logOmega2 <- c(log(omega2)[1]-diffO2, log(omega2))
omega2 <- c(0, omega2)
nbOmega2 <- nbOmega2 +1

#omega2 regular on log scale


###
### label and positions for y-axis image plot
###
eps <- -2*10^(-15)
nb0.5 <- which(omega2-0.5 >= eps)[1]
nb1 <- which(omega2-1 >= eps)[1]
nb4 <- which(omega2-1.5 >= eps)[1]
nb8 <- which(omega2-2 >= eps)[1]
myscale <- c(0,0.5,1,1.5,2)
positions <- c(logOmega2[1], logOmega2[nb0.5], logOmega2[nb1], logOmega2[nb4], logOmega2[nb8])



########### ########### ########### ###########

res1 <- NULL
for(i in phi)
{
  print(i)
  for(j in omega2)
  {
    print(j)
    res1 <- c(res1, mclapply(1:nbSimu, FUN = one.simu,
                             N = 5000,
                             sdEta = sqrt(j),
                             sdNu = 1,
                             phi = i,
                             type = "none",
                             nbSeg = 1,
                             jumpSize = 0,
                             nbK = nbK,
                             varType = "MAD",
                             mc.cores = cores)) ## mc.cores = 8
  }
}





res2 <- NULL
for(i in phi)
{
  print(i)
  for(j in omega2)
  {
    print(j)
    res2 <- c(res2, mclapply(1:nbSimu, FUN = one.simu,
                             N = 5000,
                             sdEta = sqrt(j),
                             sdNu = 1,
                             phi = i,
                             type = "rand1",
                             nbSeg = 20,
                             jumpSize = 10,
                             nbK = nbK,
                             varType = "MAD",
                             mc.cores = cores)) ## mc.cores = 8
  }
}


res3 <- NULL
for(i in phi)
{
  print(i)
  for(j in omega2)
  {
    print(j)
    res3 <- c(res3, mclapply(1:nbSimu, FUN = one.simu,
                             N = 5000,
                             sdEta = sqrt(j),
                             sdNu = 1,
                             phi = i,
                             type = "rand1",
                             nbSeg = 40,
                             jumpSize = 10,
                             nbK = nbK,
                             varType = "MAD",
                             mc.cores = cores)) ## mc.cores = 8
  }
}

########### ########### ########### ###########

df1 <- do.call(rbind, res1)
df2 <- do.call(rbind, res2)
df3 <- do.call(rbind, res3)

save(df1, file="df1.RData")
save(df2, file="df2.RData")
save(df3, file="df3.RData")


dfmean_1 <- stats::aggregate(df1, list(rep(1:(nrow(df1)%/%nbSimu+1), each = nbSimu, len = nrow(df1))), base::mean)[-1]
z1_1 <- matrix(dfmean_1$`sdEtaEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
z2_1 <- matrix(dfmean_1$`sdNuEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
z3_1 <- matrix(dfmean_1$phiEst, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)

dfsd_1 <- stats::aggregate(df1, list(rep(1:(nrow(df1)%/%nbSimu+1), each = nbSimu, len = nrow(df1))), stats::sd)[-1]
w1_1 <- matrix(dfsd_1$`sdEtaEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
w2_1 <- matrix(dfsd_1$`sdNuEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
w3_1 <- matrix(dfsd_1$phiEst, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)


dfmean_2 <- stats::aggregate(df2, list(rep(1:(nrow(df2)%/%nbSimu+1), each = nbSimu, len = nrow(df2))), base::mean)[-1]
z1_2 <- matrix(dfmean_2$`sdEtaEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
z2_2 <- matrix(dfmean_2$`sdNuEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
z3_2 <- matrix(dfmean_2$phiEst, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)

dfsd_2 <- stats::aggregate(df2, list(rep(1:(nrow(df2)%/%nbSimu+1), each = nbSimu, len = nrow(df2))), stats::sd)[-1]
w1_2 <- matrix(dfsd_2$`sdEtaEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
w2_2 <- matrix(dfsd_2$`sdNuEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
w3_2 <- matrix(dfsd_2$phiEst, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)


dfmean_3 <- stats::aggregate(df3, list(rep(1:(nrow(df3)%/%nbSimu+1), each = nbSimu, len = nrow(df3))), base::mean)[-1]
z1_3 <- matrix(dfmean_3$`sdEtaEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
z2_3 <- matrix(dfmean_3$`sdNuEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
z3_3 <- matrix(dfmean_3$phiEst, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)

dfsd_3 <- stats::aggregate(df3, list(rep(1:(nrow(df3)%/%nbSimu+1), each = nbSimu, len = nrow(df3))), stats::sd)[-1]
w1_3 <- matrix(dfsd_3$`sdEtaEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
w2_3 <- matrix(dfsd_3$`sdNuEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
w3_3 <- matrix(dfsd_3$phiEst, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)



z1_1[,1] <- 0 #do not consider evaluation for sd_eta = 0 (because of division by 0)
w1_1[,1] <- 0 #do not consider evaluation for sd_eta = 0 (because of division by 0)
z1_2[,1] <- 0 #do not consider evaluation for sd_eta = 0 (because of division by 0)
w1_2[,1] <- 0 #do not consider evaluation for sd_eta = 0 (because of division by 0)
z1_3[,1] <- 0 #do not consider evaluation for sd_eta = 0 (because of division by 0)
w1_3[,1] <- 0 #do not consider evaluation for sd_eta = 0 (because of division by 0)


