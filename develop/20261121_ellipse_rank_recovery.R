# パラメタリカバリ: 真の次元がわかっている人工データで，周辺尤度が真の次元で山になるかを確かめる ----------
# 実データで周辺尤度が次元とともに一定の幅で増え続けた(4次元まで山がない)ので，作り物かどうかを見る。
# 真のモデルは円モデル(N = 10・大きさ nu ~ N(0, 0.4^2)・誤差 sigma = 0.5)，真の次元は 2 と 3，当てる次元は 1〜4

pacman::p_load(tidyverse, cmdstanr, posterior, bridgesampling)
source("develop/20261121_ellipse_rank_datasets.R")   # prep_data()・mod

# 円モデルから非類似度を作る ------------------------------------------------------------
simulate_circle <- function(N, P, seed, sd_nu = 0.4, sigma = 0.5) {
  set.seed(seed)
  Z <- matrix(rnorm(N * P), N, P)
  Z <- sweep(Z, 2, colMeans(Z))
  Z <- Z * sqrt(N / sum(Z^2))
  nu <- rnorm(N, 0, sd_nu)
  r <- exp(nu)
  D <- as.matrix(dist(Z))
  M <- D - outer(r, rep(1, N)) + outer(rep(1, N), r)   # mu_jk = delta_jk - r_j + r_k
  d <- M + matrix(rnorm(N * N, 0, sigma), N, N)
  diag(d) <- NA
  dimnames(d) <- list(paste0("o", 1:N), paste0("o", 1:N))
  list(d = d, Z = Z, nu = nu, D = D)
}

fit_rec <- function(sim, P, label) {
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
  nu_hat <- colMeans(unclass(as.matrix(f$draws("c", format = "matrix"))))
  sig <- as.vector(f$draws("sig"))
  f$init_model_methods()
  U <- unclass(as.matrix(f$unconstrain_draws(format = "draws_matrix")))
  colnames(U) <- paste0("u", seq_len(ncol(U)))
  lb <- setNames(rep(-Inf, ncol(U)), colnames(U))
  ub <- setNames(rep(Inf, ncol(U)), colnames(U))
  lp <- function(x, data) f$log_prob(unconstrained_variables = x, jacobian = TRUE)
  b <- bridge_sampler(samples = U, log_posterior = lp, data = NULL, lb = lb, ub = ub,
                      method = "warp3", repetitions = 5, silent = TRUE)
  tibble(data = label, P = P, logml = mean(b$logml), logml_sd = sd(b$logml),
         logml_corr = mean(b$logml) + P * log(2),   # 鏡映の補正(円モデルは回転が自由なので軸の入れ替えの補正は要らない)
         divergent = sum(f$diagnostic_summary(quiet = TRUE)$num_divergent),
         sigma = mean(sig), ktau_lo = quantile(ktau, 0.025), ktau_hi = quantile(ktau, 0.975),
         ktau = mean(ktau),
         cor_dist = cor(Dhat, c(sim$D[lower.tri(sim$D)])),
         cor_nu = cor(nu_hat, sim$nu))
}

sims <- expand_grid(trueP = 2:3, rep = 1:2) %>%
  mutate(sim = map2(trueP, rep, ~ simulate_circle(10, .x, seed = 100 * .x + .y)),
         label = sprintf("真の次元%d・データ%d", trueP, rep))
res <- sims %>%
  mutate(out = map2(sim, label, function(s, l) map_dfr(1:4, ~ fit_rec(s, .x, l)))) %>%
  select(trueP, rep, out) %>% unnest(out)
saveRDS(res, "develop/20261121_ellipse_rank_recovery.rds")
