suppressMessages(pkgload::load_all(".", quiet = TRUE))
set.seed(1)
n <- 6
# 真のランク1(複素1次元)の H から作った非対称行列
z <- complex(modulus = runif(n, 1, 3), argument = runif(n, 0, pi/2))
H0 <- outer(z, Conj(z))
X0 <- Re(H0) + Im(H0)
dimnames(X0) <- list(LETTERS[1:n], LETTERS[1:n])
X <- X0; diag(X) <- NA; X[1, 2] <- X[2, 1] <- NA

cat("== 疑い2: p = n でのふるまい ==\n")
for (p in c(1, 2, n - 1, n)) {
  r <- .impute_miners(X, p = p)
  err <- max(abs(r$data[is.na(X)] - X0[is.na(X)]))
  cat(sprintf("p=%d niter=%d converged=%s 最大誤差=%.4f\n", p, r$niter, r$converged, err))
}
# 初期値そのものと p=n の結果が同じか
init <- X; for (k in which(is.na(X))) { j <- row(X)[k]; l <- col(X)[k]; init[k] <- mean(c(X[j,], X[,l]), na.rm=TRUE) }
r <- .impute_miners(X, p = n)
cat("p=n の結果と初期値の差:", max(abs(r$data[is.na(X)] - init[is.na(X)])), "\n")

cat("\n== 疑い1: ゼロを含むデータ ==\n")
# 頻度データ風: 一部をゼロにする(真値はゼロではない構造を壊さないよう, 別例)
Y <- matrix(c(NA, 5, 0, 0, 1,
              4, NA, 0, 0, 2,
              0, 0, NA, 6, 3,
              0, 1, 7, NA, 0,
              2, 1, 2, 0, NA), 5, byrow = TRUE, dimnames = list(letters[1:5], letters[1:5]))
mean_all <- function(M, j, l) mean(c(M[j,], M[,l]), na.rm = TRUE)
mean_nz  <- function(M, j, l) { v <- c(M[j,], M[,l]); v <- v[!is.na(v) & v != 0]; mean(v) }
cat("対角の初期値(現行=ゼロ込み):", round(sapply(1:5, function(j) mean_all(Y, j, j)), 3), "\n")
cat("対角の初期値(原典=非零のみ):", round(sapply(1:5, function(j) mean_nz(Y, j, j)), 3), "\n")
# 収束後に差が残るか: 現行と, 非零平均で始めた版を比べる
miners_from <- function(X, W, p, tol = 1e-4, maxit = 1000) {
  cells <- is.na(X)
  for (it in seq_len(maxit)) {
    eg <- .hfm_eigen(.hermitian(W), rotate = "none")
    U <- eg$vectors[, 1:p, drop = FALSE]
    R <- U %*% diag(eg$values[1:p], p) %*% Conj(t(U)); R <- Re(R) + Im(R)
    ch <- max(abs(R[cells] - W[cells])); W[cells] <- R[cells]
    if (ch < tol) break
  }
  list(W = W, niter = it)
}
for (p in 1:4) {
  a <- .impute_miners(Y, p = p)
  W <- Y; for (k in which(is.na(Y))) W[k] <- mean_nz(Y, row(Y)[k], col(Y)[k])
  b <- miners_from(Y, W, p)
  cat(sprintf("p=%d 現行 niter=%d / 非零版 niter=%d / 補完値の差の最大=%.5f\n",
              p, a$niter, b$niter, max(abs(a$data[is.na(Y)] - b$W[is.na(Y)]))))
}
