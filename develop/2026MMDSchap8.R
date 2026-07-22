rm(list = ls())
pacman::p_load(tidyverse)

# (8.25)以下 ----------------------------------------------------------------

Delta <- matrix(c(
  0, 5, 3, 4,
  5, 0, 2, 2,
  3, 2, 0, 1,
  4, 2, 1, 0
), ncol = 4)
W <- matrix(1, ncol = 4, nrow = 4)
Z <- matrix(c(
  -.266, .451, .016, -.200,
  -.539, .252, -.238, .524
), ncol = 2)

n <- ncol(Delta)

BZ <- function(Delta, X, W) {
  n <- NCOL(Delta)
  B_Z <- matrix(0, ncol = n, nrow = n)
  Z <- as.matrix(dist(X))
  for (i in 1:n) {
    for (j in i:n) {
      B_Z[i, j] <- -1 * W[i, j] * Delta[i, j] / Z[i, j]
      B_Z[j, i] <- B_Z[i, j]
    }
    B_Z[i, i] <- 0
    B_Z[i, i] <- -1 * sum(B_Z[i, ])
  }
  return(B_Z)
}

BZ(Delta, Z, W)

# (8.27) ------------------------------------------------------------------


V <- matrix(0, ncol = n, nrow = n)
for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    A <- matrix(0, ncol = n, nrow = n)
    A[i, i] <- A[j, j] <- 1
    A[i, j] <- A[j, i] <- -1
    V <- V + W[i, j] * A
  }
}



sigma_function <- function(D, X, W, V) {
  n <- ncol(D)
  eta2_Delta <- sum(D[upper.tri(D)]^2)
  eta2_X <- sum(diag(t(X) %*% V %*% X))
  rho <- BZ(D, X, W)
  rho_val <- sum(diag(t(X) %*% rho %*% X))

  sigma_r <- eta2_Delta + eta2_X - 2 * rho_val

  return(sigma_r)
}

eps <- 1e-6
FLG <- TRUE
max_iter <- 100
iter <- 1
sigma_old <- sigma_function(D = Delta, X = Z, W = W, V = V)
while (FLG) {
  X_u <- 1 / n * BZ(Delta, Z, W) %*% Z
  sigma_new <- sigma_function(D = Delta, X = X_u, W = W, V = V)
  print(paste(iter, sigma_new, abs(sigma_old - sigma_new)))
  if (abs(sigma_old - sigma_new) < eps) {
    FLG <- FALSE
  } else {
    sigma_old <- sigma_new
    iter <- iter + 1
  }
  if (iter > max_iter) {
    FLG <- FALSE
  }
  Z <- X_u
}
