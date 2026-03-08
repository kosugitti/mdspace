### Modern Multidimensional Scaling Chap8,Sec6
### 8.6 Majorizing Stress
rm(list = ls())
library(tidyverse)

# functions ---------------------------------------------------------------


stress_function <- function(W, Delta, Dist) {
  # W is weight matrix
  # Delta is disparity matrix
  # Dist is given Eucredian distance matrix
  ## Each matrix size should be N by N.
  N <- NROW(W)
  tmp <- 0
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      tmp <- tmp + W[i, j] * (Delta[i, j] - Dist[i, j])^2
    }
  }
  return(tmp)
}

make_B <- function(W, Delta, Dist) {
  N <- NROW(W)
  tmp <- matrix(0, nrow = N, ncol = N)
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      tmp[i, j] <- tmp[j, i] <- -W[i, j] * Dist[i, j] / Delta[i, j]
    }
  }
  diag(tmp) <- colSums(tmp) * -1
  return(tmp)
}


double_centering <- function(D){
  N <- NROW(D)
  I <- diag(1,N)
  One <- rep(1,N)
  J <- I - (1/N) * One %*% t(One)
  tmp <- -0.5 * J %*% D %*% D %*% J
  return(tmp)
}

# test --------------------------------------------------------------------

W <- matrix(1, ncol = 4, nrow = 4)
Dist <- matrix(c(0, 5, 3, 4, 5, 0, 2, 2, 3, 2, 0, 1, 4, 2, 1, 0), ncol = 4)
Z <- matrix(c(-.266, .451, .016, -.200, -.539, .252, -.238, .524), ncol = 2)
Delta <- dist(Z) %>%
  as.matrix() %>%
  round(3)

# SMACOF ------------------------------------------------------------------

k <- 0
max_iter <- 200
eps <- 1e-6
FLG <- TRUE
N <- NROW(W)
stress <- stress_function(W, Delta, Dist)

while (FLG) {
  old_stress <- stress
  k <- k + 1
  B <- make_B(W, Delta, Dist)
  X <- (1 / N) * B %*% Z
  new_Delta <- dist(X) %>%
    as.matrix() %>%
    round(3)
  stress <- stress_function(W, new_Delta, Dist)
  diff <- abs(old_stress - stress)
  if (k > max_iter | diff < eps) {
    FLG <- FALSE
  } else {
    print(paste(k, stress, diff))
    Z <- X
    Delta <- new_Delta
  }
}

X
