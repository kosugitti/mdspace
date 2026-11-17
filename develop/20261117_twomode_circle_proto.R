# 二相三元の円モデル(Okada & Imaizumi, 1997)の試作 -------------------------------------
# 式は develop/設計メモ_二相三元円モデル_20261117.md。nASYMS の出力(表1.14)と照合する
# m_jki = w_i d_jk - v_jki (r_j - r_k)，v_jki = w_i d_jk / sqrt(sum_t ((x_jt - x_kt) / u_it)^2)
# S = sqrt(mean_i S_i^2)，S_i はソースごとのストレス式2
# 方針: w_i は母数にしない(S_i は m の定数倍で変わらないので w_i = 1 として推定し，
#       最後に原典 3.7 の正規化で w_i を決める)。X・r・U を BFGS で同時に最小化する
pacman::p_load(tidyverse, pkgload)

pkgload::load_all(".")

# データ ----------------------------------------------------------------------------
raw <- utils::read.csv("tests/testthat/fixtures/okada_imaizumi_2024_table1_14.csv",
                       check.names = FALSE)
labs <- colnames(raw)[-(1:2)]
n <- length(labs)
N <- max(raw$source)
S_list <- lapply(seq_len(N), function(i) {
  m <- as.matrix(raw[raw$source == i, -(1:2)])
  dimnames(m) <- list(labs, labs)
  m
})

# セルの添字(ソースごと) -------------------------------------------------------------
make_cells <- function(S_list) {
  lapply(S_list, function(S) {
    use <- row(S) != col(S) & !is.na(S)
    list(s = S[use], j = row(S)[use], k = col(S)[use])
  })
}
cells <- make_cells(S_list)

# 母数の詰め替え: X (n x p)，r (n)，U (N x p) ------------------------------------------
unpack <- function(par, n, p, N) {
  list(X = matrix(par[seq_len(n * p)], n, p),
       r = par[n * p + seq_len(n)],
       U = matrix(par[n * p + n + seq_len(N * p)], N, p))
}

# ソース i の m(w_i = 1) ---------------------------------------------------------------
m_source <- function(X, r, u, cl) {
  dif <- X[cl$j, , drop = FALSE] - X[cl$k, , drop = FALSE]
  d <- sqrt(rowSums(dif^2))
  q <- sqrt(rowSums(sweep(dif, 2, u, "/")^2))
  v <- d / q
  list(m = d - v * (r[cl$j] - r[cl$k]), dif = dif, d = d, q = q, v = v)
}

# S(原典 3.1・3.2) ---------------------------------------------------------------------
stress_tm <- function(par, n, p, N, cells) {
  P <- unpack(par, n, p, N)
  Si2 <- vapply(seq_len(N), function(i) {
    o <- m_source(P$X, P$r, P$U[i, ], cells[[i]])
    mh <- mdspace:::.monotone(o$m, cells[[i]]$s, TRUE)
    sum((o$m - mh)^2) / sum((o$m - mean(o$m))^2)
  }, numeric(1))
  sqrt(mean(Si2))
}

# 勾配(m の式を直接微分。原典 3.10〜3.19 と同じ量を v = d / q で書き直した) -------------------
grad_tm <- function(par, n, p, N, cells) {
  P <- unpack(par, n, p, N)
  GX <- matrix(0, n, p); Gr <- numeric(n); GU <- matrix(0, N, p)
  Si2 <- numeric(N); parts <- vector("list", N)
  for (i in seq_len(N)) {
    cl <- cells[[i]]
    o <- m_source(P$X, P$r, P$U[i, ], cl)
    m <- o$m
    mh <- mdspace:::.monotone(m, cl$s, TRUE)
    A <- sum((m - mh)^2); B <- sum((m - mean(m))^2)
    Si2[i] <- A / B
    # d(S_i^2)/dm
    parts[[i]] <- list(o = o, g = 2 * ((m - mh) - (A / B) * (m - mean(m))) / B)
  }
  S <- sqrt(mean(Si2))
  idx <- factor(c(cells[[1]]$j, cells[[1]]$k), levels = seq_len(n))
  for (i in seq_len(N)) {
    cl <- cells[[i]]; o <- parts[[i]]$o
    # dS/dm = (1 / (2 N S)) d(S_i^2)/dm
    g <- parts[[i]]$g / (2 * N * S)
    u <- P$U[i, ]
    dr <- P$r[cl$j] - P$r[cl$k]
    idx <- factor(c(cl$j, cl$k), levels = seq_len(n))
    for (t in seq_len(p)) {
      a <- o$dif[, t]
      # dd/dx_jt = a/d，dq/dx_jt = a/(u_t^2 q)，dv = (dd q - d dq) / q^2
      dd <- a / o$d
      dq <- a / (u[t]^2 * o$q)
      dv <- (dd * o$q - o$d * dq) / o$q^2
      dm <- dd - dv * dr
      GX[, t] <- GX[, t] + as.vector(tapply(c(g * dm, -g * dm), idx, sum, default = 0))
      # dq/du_t = -a^2 / (u_t^3 q)，dv/du_t = -d dq / q^2
      dvu <- o$d * (a^2 / (u[t]^3 * o$q)) / o$q^2
      GU[i, t] <- sum(g * (-dvu * dr))
    }
    Gr <- Gr + as.vector(tapply(c(-g * o$v, g * o$v), idx, sum, default = 0))
  }
  c(GX, Gr, GU)
}

# 勾配の数値チェック ----------------------------------------------------------------------
p <- 2
set.seed(1)
par0 <- c(rnorm(n * p), runif(n, 0, 0.5), runif(N * p, 0.5, 1.5))
h <- 1e-6
num <- vapply(seq_along(par0), function(a) {
  e <- replace(numeric(length(par0)), a, h)
  (stress_tm(par0 + e, n, p, N, cells) - stress_tm(par0 - e, n, p, N, cells)) / (2 * h)
}, numeric(1))
grad_check <- max(abs(num - grad_tm(par0, n, p, N, cells)))
grad_check

# nASYMS の解で S を計算して照合 ------------------------------------------------------------
cf <- utils::read.csv("tests/testthat/fixtures/nasyms_twomode_nations_conf.csv")
wt <- utils::read.csv("tests/testthat/fixtures/nasyms_twomode_nations_weights.csv")
nasyms_S <- purrr::map_dfr(2:5, function(q) {
  d <- cf[cf$ndim == q, ]
  ww <- wt[wt$ndim == q, ]
  par <- c(as.matrix(d[, paste0("D", seq_len(q))]), d$radius,
           as.matrix(ww[, paste0("u", seq_len(q))]))
  tibble::tibble(ndim = q, S_recomputed = stress_tm(par, n, q, N, cells), S_nasyms = ww$S[1])
})
nasyms_S

# 最小化(初期値: 平均の対称化行列に加算定数つき Torgerson，r = 0.1，U = 1) ---------------------------
fit_tm <- function(p, par_init = NULL, maxit = 2000) {
  if (is.null(par_init)) {
    delta <- Reduce(`+`, lapply(S_list, function(S) { D <- -S; diag(D) <- NA; D - min(D, na.rm = TRUE) })) / N
    diag(delta) <- 0
    sym <- (delta + t(delta)) / 2
    X0 <- stats::cmdscale(sym, k = p, add = TRUE)$points
    par_init <- c(X0, rep(0.1, n), rep(1, N * p))
  }
  # optim の引数 par に p が部分一致しないよう，関数で包んで渡す
  fn <- function(par) stress_tm(par, n, p, N, cells)
  gr <- function(par) grad_tm(par, n, p, N, cells)
  res <- stats::optim(par_init, fn, gr, method = "BFGS",
                      control = list(maxit = maxit, reltol = 1e-12))
  list(S = res$value, par = res$par, conv = res$convergence)
}
fits <- purrr::map(2:5, fit_tm)
tibble::tibble(ndim = 2:5,
               S_proto = purrr::map_dbl(fits, "S"),
               S_nasyms = c(0.540, 0.473, 0.431, 0.400),
               conv = purrr::map_int(fits, "conv"))

# 制約つき版(方針 a) -----------------------------------------------------------------------
# u_it = exp(b_it - mean(b))：全体の幾何平均を 1 に固定(U と r の大きさの不定性を止める)
# b_it を [-L, L] に収め，同じソースの縦横比 u_it / u_is を exp(2L) = kappa 以下にする
# (楕円が線分につぶれる退化を止める)。L-BFGS-B で解く
fit_tm_con <- function(p, kappa = 10, maxit = 5000) {
  L <- log(kappa) / 2
  to_par <- function(theta) {
    b <- theta[n * p + n + seq_len(N * p)]
    c(theta[seq_len(n * p + n)], exp(b - mean(b)))
  }
  fn <- function(theta) stress_tm(to_par(theta), n, p, N, cells)
  gr <- function(theta) {
    par <- to_par(theta)
    G <- grad_tm(par, n, p, N, cells)
    iu <- n * p + n + seq_len(N * p)
    gu <- G[iu] * par[iu]
    G[iu] <- gu - mean(gu)
    G
  }
  delta <- Reduce(`+`, lapply(S_list, function(S) { D <- -S; diag(D) <- NA; D - min(D, na.rm = TRUE) })) / N
  diag(delta) <- 0
  sym <- (delta + t(delta)) / 2
  X0 <- stats::cmdscale(sym, k = p, add = TRUE)$points
  theta0 <- c(X0, rep(0.1, n), rep(0, N * p))
  lower <- c(rep(-Inf, n * p + n), rep(-L, N * p))
  upper <- c(rep(Inf, n * p + n), rep(L, N * p))
  res <- stats::optim(theta0, fn, gr, method = "L-BFGS-B", lower = lower, upper = upper,
                      control = list(maxit = maxit, factr = 1e2, pgtol = 0))
  par <- to_par(res$par)
  P <- unpack(par, n, p, N)
  # 上限に張り付いた b の数
  b <- res$par[n * p + n + seq_len(N * p)]
  at_bound <- sum(abs(abs(b) - L) < 1e-6)
  # 表示用の正規化: 布置を重心 0・二乗和 n に，半径も同じ倍率，最小の半径を 0 に
  Xc <- sweep(P$X, 2, colMeans(P$X))
  a <- sqrt(n / sum(Xc^2))
  list(S = res$value, conv = res$convergence, X = Xc * a,
       r = (P$r - min(P$r)) * a, U = P$U / a, at_bound = at_bound,
       aspect = apply(P$U, 1, function(u) max(u) / min(u)))
}
