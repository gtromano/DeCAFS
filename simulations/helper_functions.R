computeBinMes <- function(real, pred) {
  x = rep(0, tail(real, 1))
  true = x
  esti = x
  true[real] = 1
  esti[pred] = 1
  
  # computing contingency table
  contTable = table(true, esti)
  TP = contTable[2, 2]
  TN = contTable[1, 1]
  FP = contTable[1, 2]
  FN = contTable[2, 1]
  
  acc = (TN + TP) / sum(contTable)
  sens = TP / (TP + FN)
  spec = TN / (FP + TN)
  return(list(accuracy = acc, sensitivity = sens, specificity = spec))
}

# return a plot of the scenario with the signal in grey
scenarioPlot = function (Z, cps, title, n = length(Z), sdEta = 1, sdNu = 1, phi = .98) {
  set.seed(44)
  V <- dataRWAR(n = length(Z), poisParam = 0, meanGap = 0, sdEta = sdEta, sdNu = sdNu, phi = phi)$y
  #V <- rnorm(length(Z), 0, 1)
  Y = Z + V
  # jumps
  dfseg <- data.frame(x1 = cps, y1 = Y[cps], y2 = Y[cps + 1])
  df = data.frame(t = 1:length(Z), Z = Z, Y = Y)
  output = ggplot(df) +
    geom_point(aes(x = t, y = Y), col = "grey") + 
    geom_segment(aes(x = x1, xend = x1, y = y1, yend = y2), data = dfseg, arrow = arrow(length = unit(.2, "cm"))) +
    ylab("Signal") + ggtitle(title)
  output
}


# computing the f1 score from real class and predicted (coded in 0 and 1)
computeF1Score = function(real, pred, tollerance = 2.5) {
  
  real = floor(real / (tollerance * 2))
  pred = floor(pred / (tollerance * 2))
  
  x = rep(0, tail(real, 1))
  true = x
  esti = x
  true[real] = 1
  esti[pred] = 1
  
  # computing contingency table
  contTable = table(true, esti)
  if (dim(contTable)[1] == 2 && dim(contTable)[2] == 1) contTable = matrix(c(0, 0, contTable), nr = 2, nc = 2)
  TP = contTable[2, 2]
  FP = contTable[1, 2]
  FN = contTable[2, 1]
  
  precision = TP / (TP + FP)
  recall = TP / (TP + FN)
  F1 = 2 * (precision * recall) / (precision + recall)
  F1[is.nan(F1)] = 0
  F1
}


mse <- function(real, pred) mean((real - pred)^2)


# fix for AR seg function
AR1seg_func <- function (y, Kmax = 15, rho = TRUE) 
{
  l = length(y)
  if (isTRUE(rho)) {
    rho = median((diff(y, lag = 2))^2)/median(diff(y)^2) - 1
    if(rho > 1) rho <- 1
  } 
  x = y[2:l] - rho * y[1:(l - 1)]
  S = Segmentor(x, model = 2, Kmax = Kmax)
  breaks = S@breaks
  for (i in 1:Kmax) {
    for (j in 1:i) breaks[i, j] = breaks[i, j] + 1
  }
  rm(i, j)
  parameters = S@parameters
  PP = function(t) {
    x = t
    l = length(x)
    i = 2
    while (l > 2 && i < l) {
      if (x[i] == x[i - 1] + 1 && x[i] != x[i + 1] - 1) {
        x = c(x[1:(i - 1)], x[(i + 1):l])
        l = l - 1
      }
      else i = i + 1
    }
    if (l > 1 && x[l - 1] == x[l] - 1) 
      x = x[1:(l - 1)]
    x
  }
  PPbreaks = matrix(0, nrow = Kmax, ncol = Kmax, dimnames = dimnames(breaks))
  PPbreaks[1, ] = breaks[1, ]
  for (i in 2:Kmax) {
    t = PP(breaks[i, 1:(i - 1)])
    PPbreaks[i, ] = c(t, l, rep(0, Kmax - length(t) - 1))
  }
  rm(i, t)
  fMa = function(t, mu) {
    M = c()
    t = c(0, t)
    for (i in 2:length(t)) {
      M = c(M, rep(mu[i - 1], t[i] - t[i - 1]))
    }
    M
  }
  sswg = function(br, param, series) {
    sum((series - fMa(br, param))^2)
  }
  sswgseg = function(seg, seri) {
    res = c()
    for (i in 1:(Kmax)) {
      res = c(res, sswg(seg@breaks[i, 1:i], seg@parameters[i, 
                                                           1:i], seri))
    }
    res
  }
  minushalflogB = function(t, u) {
    t = t[t != 0]
    l = length(t)
    b = log(t[1]/u)
    if (l > 1) {
      for (i in 2:l) {
        b = b + log(t[i] - t[i - 1])
      }
    }
    b = -b/2
  }
  ZS = function(seg, seri) {
    u = length(seg@data)
    Kmax = seg@Kmax
    f = function(t) minushalflogB(t, u)
    wg = sswgseg(seg, seri)
    criterion = -(((u + 1):(u - Kmax + 2))/2) * log(wg) + 
      lgamma(((u + 1):(u - Kmax + 2))/2) - (0:(Kmax - 1)) * 
      log(u) + apply(seg@breaks, 1, f)
    selected = which.max(criterion)
    selected
  }
  selected = ZS(S, x)
  SelectedBreaks = breaks[selected, 1:selected]
  PPSelectedBreaks = PPbreaks[selected, ]
  PPSelectedBreaks = PPSelectedBreaks[PPSelectedBreaks != 0]
  PPselected = length(PPSelectedBreaks)
  vec1 = c(1, PPSelectedBreaks[1:(PPselected - 1)] + 1)
  vec2 = PPSelectedBreaks[1:(PPselected)]
  m = c()
  for (i in 1:PPselected) {
    m[i] = mean(y[vec1[i]:vec2[i]])
  }
  list(data = y, rho = rho, decorrelated = x, breaks = breaks, 
       PPbreaks = PPbreaks, selected = selected, SelectedBreaks = SelectedBreaks, 
       PPSelectedBreaks = PPSelectedBreaks, PPselected = PPselected, 
       PPmean = m)
}

# l2 threshold method
l2Threshold <- function(data, beta, lambda)
{
  n <- length(data)
  test <- lambda * (data[-1] - data[-n])^2
  z <- (test > beta) * (data[-1] - data[-n])
  
  chgpt <- which(z != 0)
  jump <- z[chgpt]
  
  return(list(changepoints = c(chgpt+1, n), jumpSize = jump))
}



# this function generates the change pattern. Please feel free to add more.   

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
    up = lapply(0:19, function (k)
      rep(k * jumpSize, N * 1 / 20)) %>% unlist,
    updown = lapply(0:19, function(k)
      rep((k %% 2) * jumpSize, N * 1 / 20)) %>% unlist,
    rand1 = unlist(sapply(1:20, function(i) rep(rand1Jump[i] * jumpSize, ceiling(N * rand1CP[i]))))[1:N] 
  )
}
