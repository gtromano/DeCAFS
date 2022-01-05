#' Guided Model Selection
#' 
#' This function aids the user in selecting an appropriate model for a given sequence of observations. 
#' The function goes an interactive visualization of different model fits for different choices of initial parameter estimators and l0 penalties (\code{beta}).
#' At the end, a call to the DeCAFS function is printed, while a DeCAFS wrapper is provvided.
#' 
#' @param data A vector of observations y
#'
#' @return A function, being a wrapper of DeCAFS with the selected parameter estimators.
#' @export
#'
#' @examples
#' \dontrun{
#' y <- dataRWAR(1000, sdEta = 1, sdNu = 4, phi = .4, nbSeg = 4, jumpSize = 20, type = "updown")$y
#' DeCAFSWrapper <- guidedModelSelection(y)
#' }
guidedModelSelection <- function(data) {
  cat("Guided model selection.\nThis function will return a wrapped call to DeCAFS.\n")
  
  if(!is.numeric(data)) stop("Please provvide a vector of observations.")
  
  
  # time dependency
  cat("\nStep 1, modelling the time dependency.\n")
  
  
  
  par(mfrow = c(3, 1))
  
  s_rwar <- DeCAFS(data, modelParam = estimateParameters(data))
  plot(s_rwar, main = "#1")
  mtext(side=3, line=.2, cex=0.7, paste("RWAR fit"))

  s_ar   <- DeCAFS(data, modelParam = estimateParameters(data, model = "AR"))
  plot(s_ar, main = "#2")
  mtext(side=3, line=.2, cex=0.7, paste("AR fit"))
  
  s_rw   <- DeCAFS(data, modelParam = estimateParameters(data, model = "RW"))
  plot(s_rw, main = "#3")
  mtext(side=3, line=.2, cex=0.7, paste("RW fit"))
  
  par(mfrow = c(1, 1))
  
  
  cat("Compare the estimated signals (in red) in the plots.\n(1) RWAR model fit\n(2) AR model fit\n(3) RW model fit\n")
  
  ans1 <- readline("Which is the most appropiate fit? (1/2/3):")
  estimateP <- switch (ans1,
    "1" = function(data) estimateParameters(data),
    "2" = function(data) estimateParameters(data, model = "AR"),
    "3" = function(data) estimateParameters(data, model = "RW")
  )
  
  
  cat("\nStep 2, picking an l0 penalty.\n")
  
  
  par(mfrow = c(3, 1))
  
  est <- estimateP(data)
  
  b_def <- DeCAFS(data, modelParam = est)
  plot(b_def, main = "#1")
  abline(v = b_def$changepoints, col = 4, lty = 2)
  mtext(side=3, line=.2, cex=0.7, paste("Default Beta Penalty"))
  
  b_inf <- DeCAFS(data, beta =  2 * (est$sdEta + est$sdNu) * log(length(data)), modelParam = est)
  plot(b_inf, main = "#2")
  abline(v = b_inf$changepoints, col = 4, lty = 2)
  mtext(side=3, line=.2, cex=0.7, paste("Inflated Beta Penalty"))

  b_loglog <- DeCAFS(data, beta =  2 * log(log(length(data))), modelParam = est)
  plot(b_loglog, main = "#2")
  abline(v = b_loglog$changepoints, col = 4, lty = 2)
  mtext(side=3, line=.2, cex=0.7, paste("Log-log Beta Penalty"))
  
    
  par(mfrow = c(1, 1))
  
  
  
  cat("Compare the estimated changepoints locations (in blue) in the plots.\n(1) Default BIC penalty\n(2) Inflated BIC penalty\n(3) log(log(n)) penalty\n")
  
  ans2 <- readline("Which is the most appropiate segmentation? (1/2/3):")
  penalty <- switch (ans2,
                       "1" = function(data, p) 2 * log(length(data)),
                       "2" = function(data, p) 2 * (p$sdEta + p$sdNu) * log(length(data)),
                       "3" = function(data, p) 2 * log(log(length(data)))
                     )
  
    
  cat("Guided model selection completed.\n\nRun DeCAFS with:\n")
  cat("p <- "); print(body(estimateP))
  cat("DeCAFS(data,\n")
  cat("  modelParam = p,\n")
  cat("  beta = "); print(body(penalty));
  cat(")\n")
  
  cat("\n\n")
  
  
  function(data) {
    p <- estimateP(data)
    DeCAFS(data, beta = penalty(data, p), modelParam = p)
  }
}