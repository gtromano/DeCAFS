#' Generate a piecewise constant signal of a given length
#'
#' @param N 
#' @param type Possible change scenarios for the jump structure
#' @param nbSeg Number of segments
#' @param jumpSize maximum magnitude of a change
#'
#' @return a sequence of N values for the piecewise constant signal
#'
#' @examples scenarioGenerator(1e3, "rand1")
#' 
scenarioGenerator <- function(N, type = c("none", "up", "updown", "rand1"), nbSeg = 20, jumpSize = 1) {
  #segment length

  type <- match.arg(type)

  if (type == "rand1") {
    set.seed(42)
    rand1CP <- rpois(nbSeg, lambda = 10)
    r1 <- pmax(round(rand1CP * N / sum(rand1CP)), 1) #normalisation (delete 0 values)
    s <- sum(r1)
    if(s > N)
    {
      while(sum(r1) > N)
      {
        p <- sample(x = nbSeg, size = 1)
        if(r1[p]> 1){r1[p] <- r1[p] - 1}
      }
    } else if(s < N) {
      for(i in 1:(N-s))
      {
        p <- sample(x = nbSeg, size = 1)
        r1[p] <- r1[p] + 1
      }
    }

    #jump intensity
    set.seed(43)
    rand1Jump <- runif(nbSeg, min = -1, max = 1)
  }

  switch(
    type,
    none = rep(0, N),
    up = unlist(lapply(0:(nbSeg-1), function (k) rep(k * jumpSize, N * 1 / nbSeg))),
    updown = unlist(lapply(0:(nbSeg-1), function(k) rep((k %% 2) * jumpSize, N * 1 / nbSeg))),
    rand1 = unlist(sapply(1:nbSeg, function(i) rep(rand1Jump[i] * jumpSize, r1[i])))
  )
}
