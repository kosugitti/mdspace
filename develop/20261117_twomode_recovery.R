# 二相三元の円モデル(楕円になる版)のパラメタリカバリ -----------------------------------------
# 真のモデル: m_jki = w_i d_jk - v_jki (r_j - r_k)，v_jki = d_jk / sqrt(sum_t ((x_jt - x_kt) / u_it)^2)
#   (w_i は m 全体の定数倍なので S に効かない。w_i = 1 で作る)
# シナリオ(ソース i の楕円の縦横比 = max(u_i) / min(u_i)):
#   ellipse = 縦横比 1〜3 の普通の楕円，flat = 縦横比 5〜8 の潰れがちな楕円，circle = 縦横比 1(真は円)
# 推定は縦横比の上限 kappa つき(方針 a)。合理的初期値＋乱数3本の最良
# 判定: 局所解 = 推定の S > 真の母数での S，縦横比の回復，上限への張り付き，布置・半径の回復
pacman::p_load(tidyverse, pkgload, parallel)

pkgload::load_all(".")

# ソース i の m -----------------------------------------------------------------------
m_src <- function(X, r, u, cl) {
  dif <- X[cl$j, , drop = FALSE] - X[cl$k, , drop = FALSE]
  d <- sqrt(rowSums(dif^2))
  q <- sqrt(rowSums(sweep(dif, 2, u, "/")^2))
  v <- d / q
  list(m = d - v * (r[cl$j] - r[cl$k]), dif = dif, d = d, q = q, v = v)
}

unpack <- function(par, n, p, N) {
  list(X = matrix(par[seq_len(n * p)], n, p), r = par[n * p + seq_len(n)],
       U = matrix(par[n * p + n + seq_len(N * p)], N, p))
}

stress_tm <- function(par, n, p, N, cells) {
  P <- unpack(par, n, p, N)
  Si2 <- vapply(seq_len(N), function(i) {
    o <- m_src(P$X, P$r, P$U[i, ], cells[[i]])
    mh <- mdspace:::.monotone(o$m, cells[[i]]$s, TRUE)
    sum((o$m - mh)^2) / sum((o$m - mean(o$m))^2)
  }, numeric(1))
  sqrt(mean(Si2))
}

grad_tm <- function(par, n, p, N, cells) {
  P <- unpack(par, n, p, N)
  GX <- matrix(0, n, p); Gr <- numeric(n); GU <- matrix(0, N, p)
  Si2 <- numeric(N); parts <- vector("list", N)
  for (i in seq_len(N)) {
    cl <- cells[[i]]
    o <- m_src(P$X, P$r, P$U[i, ], cl)
    m <- o$m
    mh <- mdspace:::.monotone(m, cl$s, TRUE)
    A <- sum((m - mh)^2); B <- sum((m - mean(m))^2)
    Si2[i] <- A / B
    parts[[i]] <- list(o = o, g = 2 * ((m - mh) - (A / B) * (m - mean(m))) / B)
  }
  S <- sqrt(mean(Si2))
  for (i in seq_len(N)) {
    cl <- cells[[i]]; o <- parts[[i]]$o
    g <- parts[[i]]$g / (2 * N * S)
    u <- P$U[i, ]
    dr <- P$r[cl$j] - P$r[cl$k]
    idx <- factor(c(cl$j, cl$k), levels = seq_len(n))
    for (t in seq_len(p)) {
      a <- o$dif[, t]
      dd <- a / o$d
      dq <- a / (u[t]^2 * o$q)
      dv <- (dd * o$q - o$d * dq) / o$q^2
      dm <- dd - dv * dr
      GX[, t] <- GX[, t] + as.vector(tapply(c(g * dm, -g * dm), idx, sum, default = 0))
      dvu <- o$d * (a^2 / (u[t]^3 * o$q)) / o$q^2
      GU[i, t] <- sum(g * (-dvu * dr))
    }
    Gr <- Gr + as.vector(tapply(c(-g * o$v, g * o$v), idx, sum, default = 0))
  }
  c(GX, Gr, GU)
}

# 推定(縦横比 <= kappa，全体の幾何平均 1) ------------------------------------------------------
fit_con <- function(S_list, p, kappa = 10, nrandom = 3, seed = 1) {
  n <- nrow(S_list[[1]]); N <- length(S_list)
  cells <- lapply(S_list, function(S) {
    use <- row(S) != col(S) & !is.na(S)
    list(s = S[use], j = row(S)[use], k = col(S)[use])
  })
  # u_it = exp(b_it - mean(b))，b_it = c_i + e_it
  #   c_i(ソースの大きさ)は自由，e_it(ソースの中の形)だけを [-L, L] に収める
  #   ＝ソースの中の縦横比だけが kappa 以下(ソース同士の大きさの比は縛らない)
  L <- log(kappa) / 2
  ic <- n * p + n + seq_len(N)
  ie <- n * p + n + N + seq_len(N * p)
  to_par <- function(th) {
    b <- th[ic] + th[ie]            # N x p を列優先で: e の (i, t) に c_i が足される
    c(th[seq_len(n * p + n)], exp(b - mean(b)))
  }
  iu <- n * p + n + seq_len(N * p)
  fn <- function(th) stress_tm(to_par(th), n, p, N, cells)
  gr <- function(th) {
    par <- to_par(th); G <- grad_tm(par, n, p, N, cells)
    gu <- G[iu] * par[iu]; gb <- gu - mean(gu)
    gc <- rowSums(matrix(gb, N, p))
    c(G[seq_len(n * p + n)], gc, gb)
  }
  delta <- Reduce(`+`, lapply(S_list, function(S) { D <- -S; diag(D) <- NA; D - min(D, na.rm = TRUE) })) / N
  diag(delta) <- 0
  X0 <- stats::cmdscale((delta + t(delta)) / 2, k = p, add = TRUE)$points
  starts <- list(c(X0, rep(0.1, n), rep(0, N), rep(0, N * p)))
  set.seed(seed)
  for (s in seq_len(nrandom)) {
    starts[[s + 1]] <- c(stats::rnorm(n * p), stats::runif(n, 0, 0.5), stats::rnorm(N, 0, 0.3),
                         stats::runif(N * p, -L, L))
  }
  lower <- c(rep(-Inf, n * p + n + N), rep(-L, N * p))
  upper <- c(rep(Inf, n * p + n + N), rep(L, N * p))
  fits <- lapply(starts, function(th0) {
    stats::optim(th0, fn, gr, method = "L-BFGS-B", lower = lower, upper = upper,
                 control = list(maxit = 5000, factr = 1e3))
  })
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  P <- unpack(to_par(best$par), n, p, N)
  b <- best$par[ie]
  list(S = best$value, X = P$X, r = P$r, U = P$U, conv = best$convergence,
       at_bound = sum(abs(abs(b - 0) - L) < 1e-5), cells = cells)
}

# 真の母数とデータ ------------------------------------------------------------------------
make_truth <- function(n, N, p, scenario, noise) {
  X <- matrix(stats::rnorm(n * p), n, p)
  D <- as.matrix(stats::dist(X))
  r <- stats::runif(n, 0, 0.3 * stats::median(D[upper.tri(D)]))
  aspect <- switch(scenario,
                   ellipse = stats::runif(N, 1, 3),
                   flat = stats::runif(N, 5, 8),
                   circle = rep(1, N))
  U <- t(vapply(seq_len(N), function(i) {
    s <- stats::runif(1, 0.7, 1.3)
    u <- c(s * sqrt(aspect[i]), s / sqrt(aspect[i]), rep(s, p - 2))
    if (stats::runif(1) < 0.5) u[1:2] <- u[2:1]   # どちらの次元が長軸かはソースごとにランダム
    u
  }, numeric(p)))
  S_list <- lapply(seq_len(N), function(i) {
    cl <- list(j = row(D)[row(D) != col(D)], k = col(D)[row(D) != col(D)])
    m <- m_src(X, r, U[i, ], cl)$m
    M <- matrix(NA_real_, n, n); M[row(M) != col(M)] <- m
    E <- matrix(stats::rnorm(n * n, 0, noise * stats::sd(m)), n, n)
    S <- -(M + E); diag(S) <- NA; S
  })
  list(X = X, r = r, U = U, S_list = S_list, aspect = aspect)
}

# 回復の指標(2次元: 次元の入れ替えと符号だけを許してそろえる＝U が次元に結びつくので回転は不可) ---------
align2 <- function(Xf, Uf, Xt) {
  cand <- list()
  for (perm in list(1:2, 2:1)) for (s1 in c(-1, 1)) for (s2 in c(-1, 1)) {
    Xp <- sweep(Xf[, perm], 2, c(s1, s2), "*")
    cand[[length(cand) + 1]] <- list(X = Xp, U = Uf[, perm, drop = FALSE],
                                     err = sum((Xp - Xt)^2))
  }
  cand[[which.min(vapply(cand, `[[`, numeric(1), "err"))]]
}

norm_sol <- function(X, r, U) {
  # 布置を重心 0・二乗和 n に(半径も同じ倍率)，U を幾何平均 1 に(半径を逆に)，最小の半径を 0 に
  Xc <- sweep(X, 2, colMeans(X)); a <- sqrt(nrow(X) / sum(Xc^2))
  g <- exp(mean(log(U)))
  list(X = Xc * a, r = (r - min(r)) * a / g, U = U / g)
}

one_run <- function(i, n, N, scenario, noise, kappa) {
  set.seed(10000 * n + 1000 * N + 100 * match(scenario, c("ellipse", "flat", "circle")) +
             10 * noise * 10 + i)
  tr <- make_truth(n, N, 2, scenario, noise)
  fit <- fit_con(tr$S_list, 2, kappa = kappa, seed = i)
  s_true <- stress_tm(c(tr$X, tr$r, tr$U), n, 2, N, fit$cells)
  T <- norm_sol(tr$X, tr$r, tr$U)
  F <- norm_sol(fit$X, fit$r, fit$U)
  al <- align2(F$X, F$U, T$X)
  asp_t <- apply(T$U, 1, function(u) max(u) / min(u))
  asp_f <- apply(al$U, 1, function(u) max(u) / min(u))
  # 長軸の向きが合っているか(縦横比 1.2 以上のソースだけ)
  long_ok <- (apply(T$U, 1, which.max) == apply(al$U, 1, which.max))[asp_t > 1.2]
  tibble::tibble(scenario = scenario, n = n, N = N, noise = noise, kappa = kappa, rep = i,
                 s_true = s_true, s_fit = fit$S, at_bound = fit$at_bound,
                 aspect_true = mean(asp_t), aspect_fit = mean(asp_f),
                 aspect_logerr = mean(abs(log(asp_f) - log(asp_t))),
                 long_axis_ok = if (length(long_ok)) mean(long_ok) else NA_real_,
                 conf_rmse = sqrt(mean((al$X - T$X)^2)),
                 # 真が円なら回転が決まらないので，回転も許した誤差も見る
                 conf_rmse_rot = {sv <- svd(t(F$X) %*% T$X); sqrt(mean((F$X %*% sv$u %*% t(sv$v) - T$X)^2))},
                 radius_cor = suppressWarnings(stats::cor(F$r, T$r)))
}

# 実行 -------------------------------------------------------------------------------
jobs <- tidyr::expand_grid(scenario = c("ellipse", "flat", "circle"), n = c(10, 20),
                           N = c(2, 4), noise = c(0, 0.1), kappa = 10, rep = 1:10)
res <- parallel::mclapply(seq_len(nrow(jobs)), function(k) {
  with(jobs[k, ], one_run(rep, n, N, scenario, noise, kappa))
}, mc.cores = max(1, parallel::detectCores() - 2)) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(local_min = s_fit > s_true + 1e-4)

dir.create("develop/sim_results", showWarnings = FALSE)
readr::write_csv(res, "develop/sim_results/20261117_twomode_recovery.csv")

summary_tab <- res %>%
  dplyr::group_by(scenario, n, N, noise) %>%
  dplyr::summarise(local_min = sum(local_min), at_bound = sum(at_bound > 0),
                   aspect_true = median(aspect_true), aspect_fit = median(aspect_fit),
                   aspect_logerr = median(aspect_logerr), long_axis_ok = mean(long_axis_ok, na.rm = TRUE),
                   conf_rmse = median(conf_rmse), conf_rmse_rot = median(conf_rmse_rot),
                   radius_cor = median(radius_cor),
                   .groups = "drop")
readr::write_csv(summary_tab, "develop/sim_results/20261117_twomode_recovery_summary.csv")
summary_tab
