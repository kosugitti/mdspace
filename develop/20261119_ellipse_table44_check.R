devtools::load_all(quiet = TRUE)
s <- as.matrix(read.csv("tests/testthat/fixtures/okada_imaizumi_2024_table1_19.csv", row.names = 1))
# 表4.4(書籍 p.118)の行順: Classical, Managers, Other, Proprietors, Sales, Clerical, Craft, Lower manual, Service, Primary
tab <- matrix(c(1.191,0.656,1.665,0.624, 1.025,0.047,1.290,0.394, 0.988,-1.020,1.165,1.154,
 -0.188,0.395,1.143,0.599, -0.311,-0.342,0.149,0.175, 0.215,-0.689,0.262,0.117,
 -0.453,0.956,1.175,0.821, -1.241,0.355,1.072,1.303, -0.703,-0.697,0.679,0.430,
 -0.523,0.339,1.091,0.628), 10, byrow = TRUE)
X <- tab[, 1:2]; R <- tab[, 3:4]
cat("rownames data:", rownames(s), "\n")
cat("colSums X:", colSums(X), " sumsq:", sum(X^2), "\n")
mfun <- function(X, R) {
  n <- nrow(X); M <- matrix(NA, n, n)
  for (j in 1:n) for (k in 1:n) if (j != k) {
    dif <- X[j, ] - X[k, ]; d <- sqrt(sum(dif^2))
    rj <- 1 / sqrt(sum((dif / R[j, ])^2)); rk <- 1 / sqrt(sum((dif / R[k, ])^2))
    M[j, k] <- d - d * rj + d * rk
  }
  M
}
M <- mfun(X, R)
use <- row(s) != col(s)
S <- function(M) { m <- M[use]; mh <- mdspace:::.monotone(m, s[use], TRUE); sqrt(sum((m - mh)^2) / sum((m - mean(m))^2)) }
cat("S (表4.4, そのまま):", S(M), "\n")
cat("S (転置で読む):", S(t(M)), "\n")

# 向きの一意性の確かめ(表4.6)を再現する ------------------------------------------------
# 半軸の布置 (r_j1, r_j2) を角度 a だけ回し，回した半軸を固定して点の布置だけ S を最小化する
use_m <- function(Xv, R) { M <- mfun(matrix(Xv, 10, 2), R); S(M) }
rot_check <- function(a_deg) {
  a <- a_deg * pi / 180
  Rr <- abs(R %*% matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2))
  best <- Inf
  # 表4.4 の点の布置と乱数の点の布置から始め，最良を取る
  inits <- c(list(c(X)), lapply(1:5, function(i) { set.seed(i); rnorm(20) }))
  for (z in inits) {
    o <- optim(z, use_m, R = Rr, method = "BFGS", control = list(maxit = 2000, reltol = 1e-10))
    best <- min(best, o$value)
  }
  best
}
tm <- system.time(res <- sapply(c(10, 15, 20, -10, -15, -20), rot_check))
cat("回転の角度(正=行ベクトル×回転行列):", c(10, 15, 20, -10, -15, -20), "\n")
cat("S:", round(res, 3), "\n")
cat("本の表4.6: 反時計回り 0.360 0.383 0.401 / 時計回り 0.410 0.541 0.564\n")
cat("時間(秒):", tm["elapsed"], "\n")
# 回さないときに，表4.4 の点の布置から X だけ最小化するとどこまで下がるか
o0 <- optim(c(X), use_m, R = R, method = "BFGS", control = list(maxit = 2000, reltol = 1e-10))
cat("回さない(X だけ再最小化):", round(o0$value, 4), "\n")

# nASYMS(楕円モデル・SVD の合理的初期布置・最大5次元・反復200)の2次元解(2026-11-19 にブラウザで取得) ----
# 表示: S = 1次元 0.496, 2次元 0.358, 3次元 0.254, 4次元 0.104(Not converged), 5次元 0.212
# 回した S: 反時計回り 10/15/20 度 = 0.375/0.399/0.442，時計回り = 0.388/0.408/0.463
nas <- matrix(c(-0.686,0.817,1.722,1.178, -1.181,0.519,0.901,1.089, -0.742,-0.019,1.629,1.899,
 0.258,0.394,1.276,0.958, 0.137,-1.004,0.01,0.01, -0.793,-1.004,0.419,0.001,
 0.654,0.647,1.327,1.346, 1.067,0.23,1.284,1.658, 0.644,-0.894,0.935,0.611,
 0.642,0.314,0.938,1.188), 10, byrow = TRUE)
cat("nASYMS の解: 座標の二乗和", sum(nas[, 1:2]^2), " S =", round(S(mfun(nas[, 1:2], nas[, 3:4])), 4), "(表示 0.358)\n")
# nASYMS の解から X と R をいっしょに最小化し直すとどこまで下がるか(R は対数で正に保つ)
both <- function(par) S(mfun(matrix(par[1:20], 10, 2), matrix(exp(par[21:40]), 10, 2)))
o1 <- optim(c(nas[, 1:2], log(nas[, 3:4])), both, method = "BFGS", control = list(maxit = 5000, reltol = 1e-12))
o2 <- optim(c(X, log(R)), both, method = "BFGS", control = list(maxit = 5000, reltol = 1e-12))
cat("nASYMS の解から再最小化:", round(o1$value, 4), " 表4.4 の解から再最小化:", round(o2$value, 4), "\n")
