#' Generate a piecewise constant signal of a given length
#'
#' @param n The length of the sequence of observations.
#' @param type Possible change scenarios for the jump structure
#' @param nbSeg Number of segments
#' @param jumpSize Maximum magnitude of a change
#'
#' @return a sequence of N values for the piecewise constant signal
#'
#' @examples scenarioGenerator(1e3, "rand1")
#' 
scenarioGenerator <- function(n, type = c("none", "up", "updown", "rand1"), nbSeg = 20, jumpSize = 1) {
  #segment length

  type <- match.arg(type)

  if (type == "rand1") {
    set.seed(42)
    rand1CP <- rpois(nbSeg, lambda = 10)
    r1 <- pmax(round(rand1CP * n / sum(rand1CP)), 1) #normalisation (delete 0 values)
    s <- sum(r1)
    if(s > n)
    {
      while(sum(r1) > n)
      {
        p <- sample(x = nbSeg, size = 1)
        if(r1[p]> 1){r1[p] <- r1[p] - 1}
      }
    } else if(s < n) {
      for(i in 1:(n-s))
      {
        p <- sample(x = nbSeg, size = 1)
        r1[p] <- r1[p] + 1
      }
    }

    #jump intensity
    set.seed(43)
    rand1Jump <- runif(nbSeg, min = 0.5, max = 1) * sample(c(-1,1), size = nbSeg, replace = TRUE) 
  }

  switch(
    type,
    none = rep(0, n),
    up = unlist(lapply(0:(nbSeg-1), function (k) rep(k * jumpSize, n * 1 / nbSeg))),
    updown = unlist(lapply(0:(nbSeg-1), function(k) rep((k %% 2) * jumpSize, n * 1 / nbSeg))),
    rand1 = unlist(sapply(1:nbSeg, function(i) rep(rand1Jump[i] * jumpSize, r1[i])))
  )
}
