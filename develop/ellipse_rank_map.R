# 楕円モデルの順位尤度版(develop/ellipse_rank_std.stan)の事後確率最大値を base R で求める(試作) ------------------
# circle_rank_map.R の楕円版。Stan 版と同じ対数事後密度(ヤコビアンなし・定数も含む)を，同じ母数について最大化する。
#   母数: 基準以外の対象の座標 yr((N-1) x P)，基準の対象の位置 v(P・他の平均 + exp(v))，
#         楕円の大きさ c_raw(N)，形 et(N x P・行ごとに和が 0)，m_c，s_c，s_e，sig，閾値 g(順序づけ)
#   tau = 正規化前の座標 tau_s の重心を引いて平方和 N にしたもの
#   log phi_jp = m_c + s_c c_raw_j + s_e et_jp，rho_jk = delta_jk / sqrt(sum_p (dif_p / phi_jp)^2)
#   mu_jk = delta_jk - rho_jk + rho_kj，cut = mean(mu) + sqrt(var(mu) + sig^2) * g
#   尤度: zeta ~ N(mu, sig^2) を積分した順序プロビット(タイは一次の扱い)
# 楕円の軸は座標軸にそろう(回転で変わる)ので，円モデル版のような回転の足場は使わない。鏡映は基準の対象で止める
# 下位モデル(model): "ellipse"(大きさと形)，"circle"(大きさだけ)，"equal_area"(形だけ)，"symmetric"(どちらもなし)。
# 使わない母数は推定しない(Stan では事前分布だけで動き，最大値は c_raw = 0・s_c = 0 などの端になる)。
# 対数事後密度にはその最大値での事前密度を定数として足し，Stan の lp__ と同じ値になるようにする
# 最大化は制約のない空間で行う(s_c・s_e・sig は対数，g は最初の値と間隔の対数，et は和が 0 の部分空間の正規直交基底の座標)。
# 和が 0 の部分空間への写像は等長なので，Stan の sum_to_zero_vector と同じく事前密度に余分な項は付かない

source("develop/circle_rank_map.R")   # .cr_log_diff_pnorm()・.cr_tab() を使う

# データを母数の計算に使う形にする -------------------------------------------------------------------
# x は非類似度の正方行列(対角と NA は使わない)。基準の対象は，対称部の Torgerson 解で軸ごとに重心から最も遠い対象
.er_prep <- function(x, P, model) {
  x <- as.matrix(x)
  diag(x) <- NA
  N <- nrow(x)
  D0 <- x
  D0[is.na(D0)] <- t(D0)[is.na(D0)]
  D0[is.na(D0)] <- mean(D0, na.rm = TRUE)
  diag(D0) <- 0
  Z0 <- stats::cmdscale((D0 + t(D0)) / 2, k = P)
  Z0 <- sweep(Z0, 2, colMeans(Z0))
  anchor <- apply(abs(Z0), 2, which.max)
  Z0 <- sweep(Z0, 2, sign(Z0[cbind(anchor, seq_len(P))]), `*`)
  use <- row(x) != col(x) & !is.na(x)
  vals <- sort(unique(x[use]))
  B <- stats::contr.helmert(P)
  B <- sweep(B, 2, sqrt(colSums(B^2)), `/`)            # P x (P-1)，列は正規直交で 1 と直交
  list(N = N, P = P, L = sum(use), row = row(x)[use], col = col(x)[use],
       K = length(vals), y = match(x[use], vals), anchor = anchor, Z0 = Z0, B = B,
       use_size = model %in% c("ellipse", "circle"), use_shape = model %in% c("ellipse", "equal_area"),
       model = model)
}

# 制約のない母数のベクトルと，制約のある母数の行き来 ---------------------------------------------------
# 並び: yr(列ごと)，v，[c_raw]，[z(N x (P-1)・列ごと)]，m_c，[log s_c]，[log s_e]，log sig，g の最初と間隔の対数
.er_index <- function(dat) {
  N <- dat$N
  P <- dat$P
  sizes <- c(yr = (N - 1) * P, v = P, c_raw = if (dat$use_size) N else 0,
             z = if (dat$use_shape) N * (P - 1) else 0, m_c = 1, s_c = if (dat$use_size) 1 else 0,
             s_e = if (dat$use_shape) 1 else 0, sig = 1, g = dat$K - 1)
  end <- cumsum(sizes)
  lapply(stats::setNames(names(sizes), names(sizes)), function(nm) {
    if (sizes[[nm]] == 0) integer(0) else seq(end[[nm]] - sizes[[nm]] + 1, end[[nm]])
  })
}
.er_unpack <- function(theta, dat) {
  N <- dat$N
  P <- dat$P
  ix <- .er_index(dat)
  h <- theta[ix$g]
  list(yr = matrix(theta[ix$yr], N - 1, P), v = theta[ix$v],
       c_raw = if (dat$use_size) theta[ix$c_raw] else rep(0, N),
       et = if (dat$use_shape) matrix(theta[ix$z], N, P - 1) %*% t(dat$B) else matrix(0, N, P),
       m_c = theta[ix$m_c], s_c = if (dat$use_size) exp(theta[ix$s_c]) else 0,
       s_e = if (dat$use_shape) exp(theta[ix$s_e]) else 0, sig = exp(theta[ix$sig]),
       g = cumsum(c(h[1], exp(h[-1]))))
}
.er_pack <- function(par, dat) {
  c(par$yr, par$v, if (dat$use_size) par$c_raw, if (dat$use_shape) par$et %*% dat$B,
    par$m_c, if (dat$use_size) log(par$s_c), if (dat$use_shape) log(par$s_e), log(par$sig),
    par$g[1], log(diff(par$g)))
}

# 正規化前の座標: 基準の対象は，他の対象の平均 + exp(v) --------------------------------------------------
.er_taus <- function(par, dat) {
  N <- dat$N
  P <- dat$P
  ts <- matrix(0, N, P)
  for (p in seq_len(P)) {
    a <- dat$anchor[p]
    ts[-a, p] <- par$yr[, p]
    ts[a, p] <- mean(par$yr[, p]) + exp(par$v[p])
  }
  ts
}

# 対数事後密度(Stan の target と同じ値・定数も含める)と，制約のある母数についての勾配 ------------------------
.er_logpost <- function(par, dat, grad = FALSE) {
  N <- dat$N
  P <- dat$P
  L <- dat$L
  K <- dat$K
  # 座標: 重心を引いて平方和 N に
  ts <- .er_taus(par, dat)
  C <- sweep(ts, 2, colMeans(ts))
  ss <- sum(C^2)
  kap <- sqrt(N / ss)
  tau <- kap * C
  # 半軸
  lc <- par$m_c + par$s_c * par$c_raw
  lphi <- lc + par$s_e * par$et
  iphi2 <- exp(-2 * lphi)                               # 1 / phi^2
  # モデルの値: mu = delta - rho_jk + rho_kj
  dif <- tau[dat$row, , drop = FALSE] - tau[dat$col, , drop = FALSE]
  dif2 <- dif^2
  delta <- sqrt(rowSums(dif2))
  tj <- rowSums(dif2 * iphi2[dat$row, , drop = FALSE])
  tk <- rowSums(dif2 * iphi2[dat$col, , drop = FALSE])
  rj <- delta / sqrt(tj)
  rk <- delta / sqrt(tk)
  mu <- delta - rj + rk
  mbar <- mean(mu)
  S <- sqrt(stats::var(mu) + par$sig^2)
  cut <- mbar + S * par$g
  # 順序プロビット
  eta <- mu / par$sig
  cs <- c(-Inf, cut / par$sig, Inf)
  up <- cs[dat$y + 1] - eta
  lo <- cs[dat$y] - eta
  lp_cell <- .cr_log_diff_pnorm(up, lo)
  # 事前分布(使わない母数は最大値の点 c_raw = 0・s_c = 0 などでの値を定数として足す)
  ld0 <- stats::dnorm(0, log = TRUE)
  lp <- sum(stats::dnorm(ts, log = TRUE)) + sum(par$v) +
    sum(stats::dnorm(par$c_raw, log = TRUE)) +
    sum(stats::dnorm(par$et, log = TRUE)) + N * 0.5 * log(2 * pi) +
    stats::dnorm(par$m_c, log = TRUE) +
    stats::dnorm(par$s_c, log = TRUE) + log(2) +
    stats::dnorm(par$s_e, log = TRUE) + log(2) +
    stats::dnorm(par$sig, log = TRUE) + log(2) +
    sum(stats::dnorm(par$g, log = TRUE)) + lgamma(K) +
    sum(lp_cell)
  if (!grad) return(lp)

  # 勾配: 尤度から mu・sig・g へ(circle_rank_map.R と同じ)
  g_up <- ifelse(is.finite(up), exp(stats::dnorm(up, log = TRUE) - lp_cell), 0)
  g_lo <- ifelse(is.finite(lo), -exp(stats::dnorm(lo, log = TRUE) - lp_cell), 0)
  G_eta <- -(g_up + g_lo)
  G_cs <- .cr_tab(dat$y + 1, g_up, K + 1) + .cr_tab(dat$y, g_lo, K + 1)
  G_cut <- G_cs[2:K] / par$sig
  sGc <- sum(G_cut)
  sGcg <- sum(G_cut * par$g)
  G_mu <- G_eta / par$sig + sGc / L + sGcg * (mu - mbar) / ((L - 1) * S)
  G_sig <- -sum(G_eta * mu) / par$sig^2 - sum(G_cut * cut) / par$sig + sGcg * par$sig / S - par$sig
  G_g <- G_cut * S - par$g
  # mu から dif と log phi へ
  #   d delta / d dif_p = dif_p / delta
  #   d rho_j / d dif_p = (dif_p / delta) / sqrt(tj) - delta tj^(-3/2) dif_p / phi_jp^2
  #   d rho_j / d log phi_jp = delta tj^(-3/2) dif_p^2 / phi_jp^2
  wj <- rj / tj                                         # delta tj^(-3/2)
  wk <- rk / tk
  Ij <- iphi2[dat$row, , drop = FALSE]
  Ik <- iphi2[dat$col, , drop = FALSE]
  dd <- dif / delta
  D_dif <- dd - (dd / sqrt(tj) - wj * dif * Ij) + (dd / sqrt(tk) - wk * dif * Ik)
  G_dif <- D_dif * G_mu
  Wj <- -(G_mu * wj) * dif2 * Ij                        # 行 j の log phi へ
  Wk <- (G_mu * wk) * dif2 * Ik                         # 列 k の log phi へ
  G_lphi <- matrix(0, N, P)
  G_tau <- matrix(0, N, P)
  for (p in seq_len(P)) {
    G_lphi[, p] <- .cr_tab(dat$row, Wj[, p], N) + .cr_tab(dat$col, Wk[, p], N)
    G_tau[, p] <- .cr_tab(dat$row, G_dif[, p], N) - .cr_tab(dat$col, G_dif[, p], N)
  }
  G_lc <- rowSums(G_lphi)
  G_et <- par$s_e * G_lphi - par$et
  G_se <- sum(G_lphi * par$et) - par$s_e
  G_craw <- par$s_c * G_lc - par$c_raw
  G_sc <- sum(G_lc * par$c_raw) - par$s_c
  G_mc <- sum(G_lc) - par$m_c
  # 座標: tau -> C -> tau_s(事前分布つき) -> yr・v
  G_C <- kap * (G_tau - C * sum(G_tau * C) / ss)
  G_ts <- sweep(G_C, 2, colMeans(G_C)) - ts
  G_yr <- matrix(0, N - 1, P)
  G_v <- numeric(P)
  for (p in seq_len(P)) {
    a <- dat$anchor[p]
    G_yr[, p] <- G_ts[-a, p] + G_ts[a, p] / (N - 1)
    G_v[p] <- G_ts[a, p] * exp(par$v[p]) + 1
  }
  list(lp = lp, yr = G_yr, v = G_v, c_raw = G_craw, et = G_et, m_c = G_mc, s_c = G_sc, s_e = G_se,
       sig = G_sig, g = G_g)
}

# 制約のない母数についての勾配(連鎖律) ---------------------------------------------------------------
.er_grad_theta <- function(theta, dat) {
  par <- .er_unpack(theta, dat)
  G <- .er_logpost(par, dat, grad = TRUE)
  rc <- rev(cumsum(rev(G$g)))                            # g_k = h_1 + sum_{m <= k} exp(h_m)
  Gh <- c(rc[1], rc[-1] * diff(par$g))
  c(G$yr, G$v, if (dat$use_size) G$c_raw, if (dat$use_shape) G$et %*% dat$B,
    G$m_c, if (dat$use_size) G$s_c * par$s_c, if (dat$use_shape) G$s_e * par$s_e, G$sig * par$sig, Gh)
}

# 初期値: 布置(基準の対象が正の向き)と，楕円は全員同じ大きさの円・sig 0.5・閾値は標準正規の分位 --------------
# (大きさ・形の散らばりの初期値 0.5 と m_c = -2 は Stan の試行(20261121_ellipse_rank_datasets.R)と同じ)
.er_init <- function(Z, dat) {
  N <- dat$N
  P <- dat$P
  Z <- sweep(Z, 2, colMeans(Z))
  Z <- Z * sqrt(N / sum(Z^2))
  yr <- matrix(0, N - 1, P)
  v <- numeric(P)
  for (p in seq_len(P)) {
    a <- dat$anchor[p]
    if (Z[a, p] < 0) Z[, p] <- -Z[, p]                   # 基準の対象が正の側に来るよう鏡映する
    yr[, p] <- Z[-a, p]
    v[p] <- log(Z[a, p] - mean(Z[-a, p]))
  }
  list(yr = yr, v = v, c_raw = rep(0, N), et = matrix(0, N, P), m_c = -2,
       s_c = if (dat$use_size) 0.5 else 0, s_e = if (dat$use_shape) 0.5 else 0, sig = 0.5,
       g = stats::qnorm(seq_len(dat$K - 1) / dat$K))
}

# 制約のない母数から返り値の部品を作る(元の並びのまま・回転しない) -------------------------------------
.er_result <- function(theta, dat, lab) {
  N <- dat$N
  P <- dat$P
  par <- .er_unpack(theta, dat)
  ts <- .er_taus(par, dat)
  C <- sweep(ts, 2, colMeans(ts))
  tau <- sqrt(N / sum(C^2)) * C
  phi <- exp(par$m_c + par$s_c * par$c_raw + par$s_e * par$et)
  dif <- tau[dat$row, , drop = FALSE] - tau[dat$col, , drop = FALSE]
  delta <- sqrt(rowSums(dif^2))
  mu <- delta - delta / sqrt(rowSums(dif^2 / phi[dat$row, , drop = FALSE]^2)) +
    delta / sqrt(rowSums(dif^2 / phi[dat$col, , drop = FALSE]^2))
  dn <- list(lab, paste0("D", seq_len(P)))
  dimnames(tau) <- dimnames(phi) <- dn
  list(conf = tau, semiaxis = phi, size = stats::setNames(exp(rowMeans(log(phi))), lab),
       sigma = par$sig, cut = mean(mu) + sqrt(stats::var(mu) + par$sig^2) * par$g, mu = mu,
       negative_mu = sum(mu < 0), par = par)
}

# 本体 --------------------------------------------------------------------------------------------
# x: 非類似度(type = "dissimilarity")か類似度の正方行列。ndim: 2 以上(1 次元は推定しない)。
# model: "ellipse"・"circle"・"equal_area"・"symmetric"。nrandom: 乱数の布置から始める本数(Torgerson 解に加えて)
# 返り値の conf は回さない(楕円の軸が座標軸)。semiaxis は半軸 phi(N x P)，size は半軸の幾何平均
ellipse_rank_map <- function(x, ndim = 2, type = c("dissimilarity", "similarity"),
                             model = c("ellipse", "circle", "equal_area", "symmetric"), nrandom = 10,
                             maxit = 5000, reltol = 1e-10, seed = 1) {
  type <- match.arg(type)
  model <- match.arg(model)
  x <- as.matrix(x)
  # 次元は 2 以上(MDS として布置を描くため・1 次元は距離 |x_j - x_k| の角で勾配法が止まる)
  if (length(ndim) != 1 || !is.numeric(ndim) || ndim != round(ndim) || ndim < 2 || ndim > nrow(x) - 1) {
    stop("`ndim` must be an integer between 2 and n - 1.")
  }
  if (type == "similarity") x <- -x                    # 順位しか使わないので符号を変えるだけ
  dat <- .er_prep(x, ndim, model)
  N <- dat$N
  P <- ndim
  lab <- rownames(x)
  if (is.null(lab)) lab <- paste0("o", seq_len(N))
  fn <- function(theta) {
    v <- -.er_logpost(.er_unpack(theta, dat), dat)
    if (is.finite(v)) v else 1e100                      # 直線探索で飛びすぎたときは大きな値を返す
  }
  gr <- function(theta) -.er_grad_theta(theta, dat)
  starts <- list(dat$Z0)
  if (nrandom > 0) {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(seed)
    for (i in seq_len(nrandom)) starts[[i + 1]] <- matrix(stats::rnorm(N * P), N, P)
  }
  fits <- lapply(starts, function(Z) {
    th0 <- .er_pack(.er_init(Z, dat), dat)
    tryCatch(stats::optim(th0, fn, gr, method = "BFGS", control = list(maxit = maxit, reltol = reltol)),
             error = function(e) NULL)
  })
  lps <- vapply(fits, function(f) if (is.null(f)) NA_real_ else -f$value, numeric(1))
  best <- fits[[which.max(lps)]]
  out <- .er_result(best$par, dat, lab)
  structure(c(out[c("conf", "semiaxis", "size", "sigma", "cut", "negative_mu")],
              list(lp = -best$value, lp_starts = lps, converged = best$convergence == 0,
                   niter = unname(best$counts["gradient"]), theta = best$par, dat = dat,
                   model = model, par = out$par)),
            class = "ellipse_rank_map")
}
