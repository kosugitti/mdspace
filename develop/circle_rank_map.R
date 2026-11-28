# 円モデルの順位尤度版(develop/circle_rank.stan)の事後確率最大値を base R で求める(試作) -----------------------
# Stan 版と同じ対数事後密度(ヤコビアンなし)を，同じ母数について最大化する。パッケージに移すことを考えて base R のみ。
#   母数: 足場の座標 Lx((N-1) x P の下三角・対角が正)，半径の差 a_raw(和が 0)，s_a，sig，閾値 g(順序づけ)
#   mu_jk = d_jk - a_j + a_k(a = s_a * a_raw)，cut = mean(mu) + sqrt(var(mu) + sig^2) * g
#   尤度: zeta ~ N(mu, sig^2) を積分した順序プロビット(タイは一次の扱い)
# 最大化は制約のない空間で行う(対角と s_a・sig は対数，g は最初の値と間隔の対数，a_raw は最後を残りの和の符号違いに)。
# ヤコビアンを足さないので，最大値の場所は Stan の optimize(jacobian = FALSE)と同じになる
# 解析的な勾配つきで stats::optim(BFGS)に渡す

# データを母数の計算に使う形にする -------------------------------------------------------------------
# x は非類似度の正方行列(対角と NA は使わない)。足場の P+1 個の点を先頭に並べ替える
.cr_prep <- function(x, P) {
  x <- as.matrix(x)
  diag(x) <- NA
  N <- nrow(x)
  # 対称部の Torgerson 解(欠測は転置の値で埋める)
  D0 <- x
  D0[is.na(D0)] <- t(D0)[is.na(D0)]
  diag(D0) <- 0
  Z0 <- stats::cmdscale((D0 + t(D0)) / 2, k = P)
  sc <- .cr_pick_scaffold(Z0)
  ord <- c(sc, setdiff(seq_len(N), sc))
  x <- x[ord, ord]
  Z0 <- Z0[ord, , drop = FALSE]
  use <- row(x) != col(x) & !is.na(x)
  vals <- sort(unique(x[use]))
  list(N = N, P = P, L = sum(use), row = row(x)[use], col = col(x)[use],
       K = length(vals), y = match(x[use], vals), ord = ord, Z0 = Z0)
}

# 足場の点を選ぶ: 重心から最も遠い点，そこからこれまでの点が張る空間から最も離れた点を順に ---------------------
.cr_pick_scaffold <- function(Z) {
  P <- ncol(Z)
  Zc <- sweep(Z, 2, colMeans(Z))
  idx <- which.max(rowSums(Zc^2))
  R <- sweep(Z, 2, Z[idx, ])
  for (p in seq_len(P)) {
    nxt <- which.max(rowSums(R^2))
    idx <- c(idx, nxt)
    u <- R[nxt, ] / sqrt(sum(R[nxt, ]^2))
    R <- R - tcrossprod(R %*% u, u)
  }
  idx
}

# 布置を足場の形(点 1 が原点・点 2〜N が下三角で対角が正)に回す -----------------------------------------
.cr_to_scaffold <- function(Z) {
  A <- sweep(Z[-1, , drop = FALSE], 2, Z[1, ])
  Lx <- t(qr.R(qr(t(A))))
  sweep(Lx, 2, sign(diag(Lx[seq_len(ncol(Z)), , drop = FALSE])), `*`)
}

# 制約のない母数のベクトルと，制約のある母数の行き来 ---------------------------------------------------
# 並び: Lx の下三角(列ごと・対角は対数)，a_raw の最初の N-1 個，log(s_a)，log(sig)，g の最初と間隔の対数
.cr_lower <- function(N, P) {
  which(row(matrix(0, N - 1, P)) >= col(matrix(0, N - 1, P)))
}
.cr_unpack <- function(theta, dat) {
  N <- dat$N
  P <- dat$P
  K <- dat$K
  low <- .cr_lower(N, P)
  nL <- length(low)
  Lx <- matrix(0, N - 1, P)
  Lx[low] <- theta[seq_len(nL)]
  dg <- cbind(seq_len(P), seq_len(P))
  Lx[dg] <- exp(Lx[dg])
  z <- theta[nL + seq_len(N - 1)]
  h <- theta[nL + N + 1 + seq_len(K - 1)]
  list(Lx = Lx, a_raw = c(z, -sum(z)), s_a = exp(theta[nL + N]), sig = exp(theta[nL + N + 1]),
       g = cumsum(c(h[1], exp(h[-1]))))
}
.cr_pack <- function(par, dat) {
  P <- dat$P
  Lx <- par$Lx
  dg <- cbind(seq_len(P), seq_len(P))
  Lx[dg] <- log(Lx[dg])
  c(Lx[.cr_lower(dat$N, P)], par$a_raw[-dat$N], log(par$s_a), log(par$sig),
    par$g[1], log(diff(par$g)))
}

# 対数事後密度(Stan の target と同じ値・定数も含める)と，制約のある母数についての勾配 ------------------------
.cr_logpost <- function(par, dat, grad = FALSE) {
  N <- dat$N
  P <- dat$P
  L <- dat$L
  K <- dat$K
  Lx <- par$Lx
  # 座標: 点 1 を原点に置き，重心を引いて平方和 N に
  X <- rbind(0, Lx)
  C <- sweep(X, 2, colMeans(X))
  ss <- sum(C^2)
  kap <- sqrt(N / ss)
  tau <- kap * C
  # モデルの値と閾値
  a <- par$s_a * par$a_raw
  dif <- tau[dat$row, , drop = FALSE] - tau[dat$col, , drop = FALSE]
  delta <- sqrt(rowSums(dif^2))
  mu <- delta - a[dat$row] + a[dat$col]
  mbar <- mean(mu)
  S <- sqrt(stats::var(mu) + par$sig^2)
  cut <- mbar + S * par$g
  # 順序プロビット: P(y = k) = Phi(c_k - eta) - Phi(c_{k-1} - eta)
  eta <- mu / par$sig
  cs <- c(-Inf, cut / par$sig, Inf)
  up <- cs[dat$y + 1] - eta
  lo <- cs[dat$y] - eta
  lp_cell <- .cr_log_diff_pnorm(up, lo)
  # 事前分布
  dgL <- Lx[cbind(seq_len(P), seq_len(P))]
  lp <- -0.5 * (N - 1) * P * log(2 * pi) - 0.5 * ss - 0.5 * P * log(N) +
    sum((P - seq_len(P)) * log(dgL)) +
    P * log(2) + 0.5 * P^2 * log(pi) - .cr_lmgamma(P, 0.5 * P) +
    sum(stats::dnorm(par$a_raw, log = TRUE)) + 0.5 * log(2 * pi) +
    stats::dnorm(par$s_a, log = TRUE) + log(2) +
    stats::dnorm(par$sig, log = TRUE) + log(2) +
    sum(stats::dnorm(par$g, log = TRUE)) + lgamma(K) +
    sum(lp_cell)
  if (!grad) return(lp)

  # 勾配: 尤度のセルの値の，上端・下端の引数についての微分
  g_up <- ifelse(is.finite(up), exp(stats::dnorm(up, log = TRUE) - lp_cell), 0)
  g_lo <- ifelse(is.finite(lo), -exp(stats::dnorm(lo, log = TRUE) - lp_cell), 0)
  G_eta <- -(g_up + g_lo)
  # cs の各要素について(端の 2 つ = -Inf・Inf は捨てる)
  G_cs <- .cr_tab(dat$y + 1, g_up, K + 1) + .cr_tab(dat$y, g_lo, K + 1)
  G_cut <- G_cs[2:K] / par$sig
  # mu について(直接の分と，閾値の平均・広がりを通る分)
  sGc <- sum(G_cut)
  sGcg <- sum(G_cut * par$g)
  G_mu <- G_eta / par$sig + sGc / L + sGcg * (mu - mbar) / ((L - 1) * S)
  # sig について(eta・cut/sig の直接の分と，S を通る分)＋事前分布
  G_sig <- -sum(G_eta * mu) / par$sig^2 - sum(G_cut * cut) / par$sig + sGcg * par$sig / S - par$sig
  G_g <- G_cut * S - par$g
  # 半径の差
  G_a <- -.cr_tab(dat$row, G_mu, N) + .cr_tab(dat$col, G_mu, N)
  G_araw <- par$s_a * G_a - par$a_raw
  G_sa <- sum(G_a * par$a_raw) - par$s_a
  # 座標: delta -> tau -> C -> X -> Lx
  u <- dif * (G_mu / delta)
  G_tau <- matrix(0, N, P)
  for (p in seq_len(P)) {
    G_tau[, p] <- .cr_tab(dat$row, u[, p], N) - .cr_tab(dat$col, u[, p], N)
  }
  G_C <- kap * (G_tau - C * sum(G_tau * C) / ss) - C
  G_X <- sweep(G_C, 2, colMeans(G_C))
  G_Lx <- G_X[-1, , drop = FALSE]
  G_Lx[cbind(seq_len(P), seq_len(P))] <- G_Lx[cbind(seq_len(P), seq_len(P))] + (P - seq_len(P)) / dgL
  list(lp = lp, Lx = G_Lx, a_raw = G_araw, s_a = G_sa, sig = G_sig, g = G_g)
}

# 制約のない母数についての勾配(連鎖律) ---------------------------------------------------------------
.cr_grad_theta <- function(theta, dat) {
  par <- .cr_unpack(theta, dat)
  G <- .cr_logpost(par, dat, grad = TRUE)
  P <- dat$P
  N <- dat$N
  GL <- G$Lx
  dg <- cbind(seq_len(P), seq_len(P))
  GL[dg] <- GL[dg] * par$Lx[dg]                         # 対角は対数
  Gz <- G$a_raw[-N] - G$a_raw[N]                       # a_raw[N] = -sum(z)
  Gg <- G$g
  rc <- rev(cumsum(rev(Gg)))                            # g_k = h_1 + sum_{m <= k} exp(h_m)
  Gh <- c(rc[1], rc[-1] * diff(par$g))
  c(GL[.cr_lower(N, P)], Gz, G$s_a * par$s_a, G$sig * par$sig, Gh)
}

# 補助: log(Phi(u) - Phi(l))(u > l)を桁落ちなく ------------------------------------------------------
.cr_log_diff_pnorm <- function(u, l) {
  out <- numeric(length(u))
  pos <- !is.na(l) & l > 0                             # 両方が右の裾なら上側確率の差で
  # (直線探索で母数が飛んで NaN になったときは NaN のまま返す＝optim の BFGS は歩幅を縮めてやり直す)
  out[pos] <- log(stats::pnorm(l[pos], lower.tail = FALSE) - stats::pnorm(u[pos], lower.tail = FALSE))
  out[!pos] <- log(stats::pnorm(u[!pos]) - stats::pnorm(l[!pos]))
  # 端のセル(片方が無限大)は 1 つの裾の対数確率で
  e1 <- is.infinite(l)
  out[e1] <- stats::pnorm(u[e1], log.p = TRUE)
  e2 <- is.infinite(u)
  out[e2] <- stats::pnorm(l[e2], lower.tail = FALSE, log.p = TRUE)
  out
}
.cr_tab <- function(idx, w, n) {
  out <- numeric(n)
  s <- rowsum(w, idx)
  out[as.integer(rownames(s))] <- s[, 1]
  out
}
.cr_lmgamma <- function(P, x) {
  0.25 * P * (P - 1) * log(pi) + sum(lgamma(x + (1 - seq_len(P)) / 2))
}

# 初期値: 布置(足場の形)と，半径の差 0・s_a 0.3・sig 0.5・閾値は標準正規の分位 ----------------------------
.cr_init <- function(Z, dat) {
  Z <- sweep(Z, 2, colMeans(Z))
  Z <- Z * sqrt(dat$N / sum(Z^2))
  list(Lx = .cr_to_scaffold(Z), a_raw = rep(0, dat$N), s_a = 0.3, sig = 0.5,
       g = stats::qnorm(seq_len(dat$K - 1) / dat$K))
}

# sigma の周辺事後分布のラプラス近似 ----------------------------------------------------------------
# 同時の最大値では sigma が小さく出る(他の母数が多いと，sigma の小さい側の細く高い山が選ばれる)。
# そこで sigma を固定した山を正規分布で近似して他の母数を積分し，log sigma の周辺事後密度
#   log p(log sig | y) ≈ log p(theta_hat, log sig | y) + log|ヤコビアン| - 0.5 log det H + 定数
# を最大にする sigma を求める。積分は制約のない空間で行うので，変数変換のヤコビアンを足す

# 制約のない母数のうち，対数を取った位置(対角・s_a・sig・閾値の間隔)。ヤコビアンの対数はこれらの和
.cr_jac_index <- function(dat) {
  N <- dat$N
  P <- dat$P
  nL <- length(.cr_lower(N, P))
  dpos <- match((seq_len(P) - 1) * (N - 1) + seq_len(P), .cr_lower(N, P))
  c(dpos, nL + N, nL + N + 1, nL + N + 1 + seq_len(dat$K - 1)[-1])
}
.cr_sig_index <- function(dat) length(.cr_lower(dat$N, dat$P)) + dat$N + 1

# sigma(対数)を固定したときの，他の母数についての負の対数密度と勾配(jac = TRUE でヤコビアンを足す)
.cr_fix_sig <- function(ls, dat, jac) {
  isig <- .cr_sig_index(dat)
  jidx <- setdiff(.cr_jac_index(dat), isig)
  full <- function(v) {
    th <- numeric(length(v) + 1)
    th[-isig] <- v
    th[isig] <- ls
    th
  }
  fn <- function(v) {
    th <- full(v)
    -(.cr_logpost(.cr_unpack(th, dat), dat) + if (jac) sum(th[jidx]) + ls else 0)
  }
  gr <- function(v) {
    g <- .cr_grad_theta(full(v), dat)
    if (jac) g[jidx] <- g[jidx] + 1
    -g[-isig]
  }
  list(fn = fn, gr = gr, full = full)
}

# theta_map から始めて，log sigma の周辺事後密度(ラプラス近似)を 1 次元で最大化する
.cr_sigma_laplace <- function(theta_map, dat, width = log(3), maxit = 5000, reltol = 1e-10) {
  isig <- .cr_sig_index(dat)
  ls0 <- theta_map[isig]
  cur <- theta_map[-isig]
  trace <- NULL
  marg <- function(ls) {
    f <- .cr_fix_sig(ls, dat, jac = TRUE)
    r <- stats::optim(cur, f$fn, f$gr, method = "BFGS", control = list(maxit = maxit, reltol = reltol))
    cur <<- r$par                                      # 次の候補はこの答えから始める
    H <- stats::optimHess(r$par, f$fn, f$gr)
    H <- (H + t(H)) / 2
    ld <- determinant(H, logarithm = TRUE)
    val <- if (ld$sign > 0) -r$value - 0.5 * as.numeric(ld$modulus) + 0.5 * length(cur) * log(2 * pi) else -Inf
    trace <<- rbind(trace, data.frame(sigma = exp(ls), log_marg = val, log_height = -r$value))
    val
  }
  lower <- ls0 - 0.1
  upper <- ls0 + width
  opt <- stats::optimize(marg, c(lower, upper), maximum = TRUE, tol = 1e-3)
  if (upper - opt$maximum < 0.01) warning("sigma の周辺事後分布の山が探した範囲の上端にあります。width を広げてください")
  trace <- trace[order(trace$sigma), ]
  rownames(trace) <- NULL
  list(log_sigma = opt$maximum, log_marg = opt$objective, trace = trace)
}

# 制約のない母数から返り値の部品を作る(元の並びに戻し，布置は主軸に回す。半径は最小を 0 に・OI と同じ) -----------
.cr_result <- function(theta, dat, lab) {
  N <- dat$N
  P <- dat$P
  par <- .cr_unpack(theta, dat)
  X <- rbind(0, par$Lx)
  C <- sweep(X, 2, colMeans(X))
  tau <- sqrt(N / sum(C^2)) * C
  conf <- matrix(NA_real_, N, P)
  conf[dat$ord, ] <- tau %*% svd(tau)$v
  a <- par$s_a * par$a_raw
  radius <- numeric(N)
  radius[dat$ord] <- a - min(a)
  rownames(conf) <- names(radius) <- lab
  mu <- sqrt(rowSums((tau[dat$row, , drop = FALSE] - tau[dat$col, , drop = FALSE])^2)) -
    a[dat$row] + a[dat$col]
  list(conf = conf, radius = radius, sigma = par$sig,
       cut = mean(mu) + sqrt(stats::var(mu) + par$sig^2) * par$g, par = par)
}

# 本体 --------------------------------------------------------------------------------------------
# x: 非類似度(type = "dissimilarity")か類似度の正方行列。ndim: 2 以上(1 次元は推定しない)。nrandom: 乱数の布置から始める本数(Torgerson 解に加えて)
# sigma = "map"(既定): 同時の最大値のまま。sigma は小さめに出る(N=12 で MCMC の約 0.8 倍，N=50 で約 0.95 倍)
# sigma = "laplace": sigma を周辺事後分布(ラプラス近似)の山で置き直し，その sigma で他の母数を最大化し直す
#   (sigma は MCMC に近づくが半径が膨らむ・20261124_circle_rank_laplace_check.R)
circle_rank_map <- function(x, ndim = 2, type = c("dissimilarity", "similarity"), nrandom = 10,
                            sigma = c("map", "laplace"), maxit = 5000, reltol = 1e-10, seed = 1) {
  type <- match.arg(type)
  sigma <- match.arg(sigma)
  x <- as.matrix(x)
  # 次元は 2 以上(MDS として布置を描くため・1 次元は距離 |x_j - x_k| の角で勾配法が止まる)。
  # 足場に P+1 個の点が要るので N - 1 以下
  if (length(ndim) != 1 || !is.numeric(ndim) || ndim != round(ndim) || ndim < 2 || ndim > nrow(x) - 1) {
    stop("`ndim` must be an integer between 2 and n - 1.")
  }
  if (type == "similarity") x <- -x                    # 順位しか使わないので符号を変えるだけ
  dat <- .cr_prep(x, ndim)
  N <- dat$N
  P <- ndim
  lab <- rownames(x)
  if (is.null(lab)) lab <- paste0("o", seq_len(N))
  fn <- function(theta) -.cr_logpost(.cr_unpack(theta, dat), dat)
  gr <- function(theta) -.cr_grad_theta(theta, dat)
  # 初期値: Torgerson 解と，乱数の布置 nrandom 本
  starts <- list(dat$Z0)
  if (nrandom > 0) {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(seed)
    for (i in seq_len(nrandom)) starts[[i + 1]] <- matrix(stats::rnorm(N * P), N, P)
  }
  fits <- lapply(starts, function(Z) {
    th0 <- .cr_pack(.cr_init(Z, dat), dat)
    tryCatch(stats::optim(th0, fn, gr, method = "BFGS", control = list(maxit = maxit, reltol = reltol)),
             error = function(e) NULL)
  })
  lps <- vapply(fits, function(f) if (is.null(f)) NA_real_ else -f$value, numeric(1))
  best <- fits[[which.max(lps)]]
  theta <- best$par
  sigma_map <- .cr_unpack(theta, dat)$sig
  lap <- NULL
  if (sigma == "laplace") {
    lap <- .cr_sigma_laplace(theta, dat, maxit = maxit, reltol = reltol)
    # 置き直した sigma で，他の母数を最大化し直す(ヤコビアンなし＝sigma を与えたときの事後確率最大値)
    f <- .cr_fix_sig(lap$log_sigma, dat, jac = FALSE)
    r <- stats::optim(theta[-.cr_sig_index(dat)], f$fn, f$gr, method = "BFGS",
                      control = list(maxit = maxit, reltol = reltol))
    theta <- f$full(r$par)
  }
  out <- .cr_result(theta, dat, lab)
  structure(c(out[c("conf", "radius", "sigma", "cut")],
              list(sigma_map = sigma_map, sigma_trace = lap$trace, log_marg_sigma = lap$log_marg,
                   lp = max(lps, na.rm = TRUE), lp_starts = lps,
                   n_same = sum(max(lps, na.rm = TRUE) - lps < 0.1, na.rm = TRUE),
                   convergence = best$convergence, niter = best$counts[["function"]],
                   par = out$par, theta = theta, dat = dat, call = match.call())),
            class = "mdspace_circle_rank")
}
