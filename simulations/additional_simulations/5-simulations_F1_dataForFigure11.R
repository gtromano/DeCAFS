rm(list=ls())

# LIBRARIES
#devtools::install_github("gtromano/DeCAFS")
#devtools::install_github("vrunge/ARRWestim")
library(DeCAFS)
library(parallel)
library(fields)


###
### NB cores for parallel computing
###
cores <- detectCores()
cores <- 60
#cores <- 8 ###CHANGE

one.simu.2stages <- function(i, N = 5*10^3, sdEta = 0.1, sdNu = 0.3, phi = 0.2,
                     type = "rand1", nbSeg = 10,
                     jumpSize = 2, nbK = 10, varType = "MAD")
{
  y <- ARRWestim::dataRWAR(N = N,
                sdEta = sdEta, sdNu = sdNu, phi = phi,
                type = type,
                nbSeg = nbSeg, jumpSize = jumpSize,
                seed = sample(1e6,1))
  #1st estimate
  res <- bestParameters(y$y, nbK = nbK, type = varType)
  #DECAFS
  deca <- DeCAFS::DeCAFS(y$y, 2*log(N),
      list(sdEta = res$EtaOpt, sdNu = res$NuOpt, phi = res$argmin))


  chpts <- c(0, deca$changepoints, N)
  K <- length(chpts) - 1
  w <- diff(chpts)

  newRes <- matrix(0, nrow = 0, ncol = 5) #start + end + 3 parameters estimated
  for(i in 2:length(chpts))
  {
    if(w[i-1] > 50){newRes <- rbind(newRes, c(chpts[(i-1):i],0,0,0))}
  }
  if(nrow(newRes) > 0)
  {
    for(i in 1:nrow(newRes))
    {
      seg <- y$y[(newRes[i,1]+1):newRes[i,2]]
      newRes[i,3:5] <- unname(unlist(bestParameters(seg, nbK = nbK, type = varType)))
    }
  }
  finalRes <- apply(newRes, weighted.mean, newRes[,2] - newRes[,1], MARGIN = 2)
  #ind, n, phi, sdEta2, sdNu2, nbK, poiP, meanGap, method, mse, beta, K
  df <- data.frame(numeric(0), numeric(0), numeric(0), numeric(0), numeric(0),
                   numeric(0), numeric(0), numeric(0), numeric(0), numeric(0),
                   numeric(0) ,stringsAsFactors = FALSE)
  colnames(df) <- c("index", "n", "sdEta", "sdNu", "phi", "nbK",
                    "nbSeg", "jumpSize", "sdEtaEst%", "sdNuEst%", "phiEst_error")
  df[1,] <- c(i,N, sdEta, sdNu, phi, nbK, nbSeg, jumpSize,
              (finalRes[3] - sdEta)/sdEta, (finalRes[4]-sdNu)/sdNu, finalRes[5] - phi)
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
    res1 <- c(res1, mclapply(1:nbSimu, FUN = one.simu.2stages,
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

df_2stage <- do.call(rbind, res1)
save(df_2stage, file="df_2stage.RData")

dfmean_1 <- stats::aggregate(df_2stage, list(rep(1:(nrow(df_2stage)%/%nbSimu+1), each = nbSimu, len = nrow(df_2stage))), base::mean)[-1]
z1 <- matrix(dfmean_1$`sdEtaEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
z2 <- matrix(dfmean_1$`sdNuEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
z3 <- matrix(dfmean_1$phiEst, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)

dfsd_1 <- stats::aggregate(df_2stage, list(rep(1:(nrow(df_2stage)%/%nbSimu+1), each = nbSimu, len = nrow(df_2stage))), stats::sd)[-1]
w1 <- matrix(dfsd_1$`sdEtaEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
w2 <- matrix(dfsd_1$`sdNuEst%`, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)
w3 <- matrix(dfsd_1$phiEst, nrow = nbPhi, ncol = nbOmega2, byrow = TRUE)

z1[,1] <- 0 #do not consider evaluation for sd_eta = 0 (because of division by 0)
w1[,1] <- 0 #do not consider evaluation for sd_eta = 0 (because of division by 0)

