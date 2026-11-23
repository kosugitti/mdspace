# パラメタリカバリ(Okada & Imaizumi 1987 の 4.1 と同じ条件) ----------------------------------
# 11/21 のリカバリ(N = 10・sd_nu 0.4・sigma 0.5)は，原典の定義で誤差 E = 0.22〜0.59・比 R = 0.36〜1.71 と
# 原典の最も厳しい条件(E = 0.3・R = 0.5)より厳しかった。原典と同じ条件にそろえ，
# 順位尤度版(円モデル)の周辺尤度・順位の一致率・sigma と，circle_mds() の stress の肘を並べる
# 原典: 対象 12・真の次元 2・R = 0.1, 0.3, 0.5・E = 0, 0.1, 0.3・各 5 回・3〜1 次元で分析
#   R = 2 n^2 var(r) / sum_{j != k} (d_jk - dbar)^2 (式 4-1 の分母の 2 項の比)
#   E = Ve / var(m)(誤差の分散と m の分散の比)，s_jk = m_jk + e_jk(E = 0 では s = 1.8 m^2 + 5.5)
#   (2026-11-22) 原典は E を分散の比と書くが，分散の比で作ると stress が表4.1 よりずっと大きい(E = 0.1 の 2 次元で 0.21 対 0.06)。
#   SD の比(誤差の SD = E * SD(m))で作ると表4.1 とほぼ一致する(E = 0.1: 0.055〜0.059 対 0.055〜0.069，E = 0.3: 0.195〜0.201 対 0.197〜0.202)ので，SD の比として使う
# ここでは E = 0.1, 0.3 だけを使う(誤差なしは順序プロビットで sigma -> 0 になり推定が壊れる)。繰り返しはまず 2 回

pacman::p_load(tidyverse, cmdstanr, posterior, bridgesampling)
source("develop/20261121_ellipse_rank_datasets.R")   # prep_data()・mod
devtools::load_all(quiet = TRUE)                     # circle_mds()

# 原典の条件で非類似度を作る ------------------------------------------------------------------
simulate_oi <- function(N = 12, P = 2, R, E, seed) {
  set.seed(seed)
  Z <- matrix(rnorm(N * P), N, P)
  Z <- sweep(Z, 2, colMeans(Z))
  Z <- Z * sqrt(N / sum(Z^2))
  D <- as.matrix(dist(Z))
  off <- row(D) != col(D)
  # 半径: 一様乱数を R に合うように伸縮し，最小を 0 に
  r0 <- runif(N)
  r0 <- r0 - mean(r0)
  denom <- sum((D[off] - mean(D[off]))^2)
  r <- r0 * sqrt(R * denom / (2 * N^2 * mean(r0^2)))   # var は母分散(sum_{j != k}(r_k - r_j)^2 = 2 n^2 var(r))
  r <- r - min(r)
  M <- D - outer(r, rep(1, N)) + outer(rep(1, N), r)   # m_jk = d_jk - r_j + r_k
  e <- matrix(rnorm(N * N, 0, E * sd(M[off])), N, N)   # E は SD の比として使う(下の注)
  d <- M + e
  diag(d) <- NA
  dimnames(d) <- list(paste0("o", 1:N), paste0("o", 1:N))
  list(d = d, Z = Z, r = r, D = D)
}

# 順位尤度版(円モデル)を当てる -----------------------------------------------------------------
fit_bayes <- function(sim, P) {
  pr <- prep_data(sim$d, "dissimilarity", P = P)
  od <- tempfile("fit_")
  dir.create(od)
  f <- mod$sample(data = modifyList(pr$dat, list(use_size = 1, use_shape = 0)),
                  init = pr$make_init, chains = 4, parallel_chains = 4, iter_warmup = 1000,
                  iter_sampling = 2500, seed = 1, refresh = 0, show_messages = FALSE,
                  output_dir = od, sig_figs = 18)
  N <- pr$dat$N
  mu <- unclass(as.matrix(f$draws("mu", format = "matrix")))
  idx <- round(seq(1, nrow(mu), length.out = 1000))
  ktau <- vapply(idx, function(i) cor(mu[i, ], pr$dat$y, method = "kendall"), numeric(1))
  tau <- unclass(as.matrix(f$draws("tau", format = "matrix")))
  Dhat <- colMeans(t(apply(tau, 1, function(v) c(dist(matrix(v, N, P))))))
  rad <- exp(colMeans(unclass(as.matrix(f$draws("c", format = "matrix")))))
  f$init_model_methods()
  U <- unclass(as.matrix(f$unconstrain_draws(format = "draws_matrix")))
  colnames(U) <- paste0("u", seq_len(ncol(U)))
  lb <- setNames(rep(-Inf, ncol(U)), colnames(U))
  ub <- setNames(rep(Inf, ncol(U)), colnames(U))
  lp <- function(x, data) f$log_prob(unconstrained_variables = x, jacobian = TRUE)
  b <- bridge_sampler(samples = U, log_posterior = lp, data = NULL, lb = lb, ub = ub,
                      method = "warp3", repetitions = 3, silent = TRUE)
  tibble(logml = mean(b$logml) + P * log(2),   # 鏡映の補正(円モデルは回転が自由なので入れ替えの補正は要らない)
         logml_sd = sd(b$logml),
         divergent = sum(f$diagnostic_summary(quiet = TRUE)$num_divergent),
         sigma = mean(f$draws("sig")),
         ktau_lo = quantile(ktau, 0.025), ktau_hi = quantile(ktau, 0.975), ktau = mean(ktau),
         cor_dist = cor(Dhat, c(sim$D[lower.tri(sim$D)])),
         cor_r = cor(rad, sim$r))
}

# circle_mds() の stress(原典と同じく 3 次元から下ろす) ------------------------------------------
fit_stress <- function(sim, P) {
  out <- circle_mds(sim$d, ndim = P, type = "dissimilarity", ndim_max = 3)
  tibble(stress = out$stress,
         cor_dist_mds = cor(c(dist(out$conf)), c(sim$D[lower.tri(sim$D)])),
         cor_r_mds = cor(out$radius, sim$r))
}

grid <- expand_grid(R = c(0.1, 0.3, 0.5), E = c(0.1, 0.3), rep = 1:2) %>%
  mutate(seed = 1000 + row_number())
res <- grid %>%
  mutate(out = pmap(list(R, E, seed), function(R, E, seed) {
    sim <- simulate_oi(R = R, E = E, seed = seed)
    map_dfr(1:3, function(P) bind_cols(tibble(P = P), fit_bayes(sim, P), fit_stress(sim, P)))
  })) %>%
  unnest(out)
saveRDS(res, "develop/20261122_ellipse_rank_recovery_oi1987.rds")
