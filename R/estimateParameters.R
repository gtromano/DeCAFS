estimateParameters = function(y) {
  N = length(y)
  y_0 = y[1:(N - 2)]
  y_1 = y[2:(N - 1)]
  y_2 = y[3:(N)]
  v1 = mad(y_1 - y_0) ^ 2
  v2 = mad(y_2 - y_0) ^ 2
  sigmaXsq = abs(v2 - v1)
  sigmaYsq = abs(v1 - v2/2)
  return(list(sigmaX = sqrt(sigmaXsq), sigmaY = sqrt(sigmaYsq)))
}

