# 円モデルの順位尤度版: sigma の事後確率最大値が MCMC の事後平均より小さい理由を確かめる -----------------------
# 20261123_circle_rank_approx.R では，sigma の事後確率最大値が MCMC の事後平均の 0.5〜0.8 倍だった。
# 考えられる理由は 2 つ:
#   (1) sigma の周辺事後分布が右に裾を引いていて，平均と最頻値がずれているだけ
#   (2) 値がすべて違うデータでは閾値が K-1 = L-1 個あり，セルの数とほぼ同じ数の母数を同時に最大化するので
#       sigma が小さく出る(母数の数がデータとともに増えるときの分散の最尤推定量と同じ仕組み・Neyman & Scott 1948)
# (1) は MCMC の周辺分布の中央値・最頻値と比べる。(2) は非類似度を分位で K 段階に粗くして閾値の数を減らし，
# 事後確率最大値と MCMC の差が縮むかを見る(順位しか使わないので，粗くするとタイが増えるだけ)

pacman::p_load(tidyverse, cmdstanr, posterior)
mod <- cmdstanr::cmdstan_model("develop/circle_rank.stan")

load_defs <- function(file, names) {
  for (e in parse(file)) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]]) %in% names) {
      eval(e, envir = globalenv())
    }
  }
}
load_defs("develop/20261122_ellipse_rank_recovery_oi1987.R", "simulate_oi")
load_defs("develop/20261123_circle_rank_scaffold.R", c("pick_scaffold", "to_scaffold", "prep_circle"))

# 非類似度を分位で K 段階に粗くする(K = Inf ならそのまま) -------------------------------------------
coarsen <- function(d, K) {
  if (is.infinite(K)) return(d)
  off <- !is.na(d)
  br <- quantile(d[off], seq(0, 1, length.out = K + 1))
  d[off] <- as.numeric(cut(d[off], br, include.lowest = TRUE))
  d
}

# MCMC と事後確率最大値の sigma を比べる -------------------------------------------------------------
fit_sigma <- function(sim, P, K) {
  d <- coarsen(sim$d, K)
  pr <- prep_circle(d, P)
  N <- pr$dat$N
  init <- function(chain_id) {
    i0 <- pr$make_init(chain_id)
    i0[c("c_raw", "m_c", "s_c")] <- NULL
    i0$a_raw <- rep(0, N)
    i0$s_a <- 0.3
    i0
  }
  od <- tempfile("fit_")
  dir.create(od)
  D_true <- sim$D[pr$ord, pr$ord][lower.tri(diag(N))]
  f <- mod$sample(data = pr$dat, init = init, chains = 4, parallel_chains = 4, iter_warmup = 1000,
                  iter_sampling = 2500, seed = 1, refresh = 0, show_messages = FALSE,
                  adapt_delta = 0.99, output_dir = od)
  sig <- c(f$draws("sig"))
  dn <- density(log(sig))                             # 最頻値は対数の目盛りの密度から(0 の近くの境界を避ける)
  tau <- unclass(as.matrix(f$draws("tau", format = "matrix")))
  D_mcmc <- colMeans(t(apply(tau, 1, function(v) c(dist(matrix(v, N, P))))))
  # 事後確率最大値: Torgerson 解の近くから 5 本の最良
  opts <- map(1:5, function(i) {
    mod$optimize(data = pr$dat, init = list(init(i)), seed = i, jacobian = FALSE, iter = 20000,
                 refresh = 0, show_messages = FALSE, output_dir = od)
  })
  lps <- vapply(opts, function(o) tryCatch(o$lp(), error = function(e) NA_real_), numeric(1))
  best <- opts[[which.max(lps)]]
  s_map <- unname(best$mle("sig"))
  D_map <- c(dist(matrix(best$mle("tau"), N, P)))
  tibble(K_used = pr$dat$K, L = pr$dat$L,
         sig_map = s_map, sig_mean = mean(sig), sig_median = median(sig),
         sig_mode = exp(dn$x[which.max(dn$y)]),       # log(sigma) の密度の最頻値
         sig_q05 = quantile(sig, 0.05), sig_q95 = quantile(sig, 0.95),
         q_of_map = mean(sig < s_map),                # 事後確率最大値が周辺分布の下から何割の所か
         ratio = s_map / mean(sig),
         cor_D = cor(D_map, D_mcmc),
         cor_D_true_map = cor(D_map, D_true), cor_D_true_mcmc = cor(D_mcmc, D_true),
         divergent = sum(f$diagnostic_summary(quiet = TRUE)$num_divergent))
}

# データと条件 --------------------------------------------------------------------------------------
dsets <- tribble(
  ~label,               ~N, ~R,  ~E,  ~seed,
  "N12・E0.5・データ1", 12, 0.3, 0.5, 2001,
  "N10・E0.7・データ2", 10, 0.3, 0.7, 2010,
  "N12・E0.3(対照)",   12, 0.3, 0.3, 1007
)
res <- expand_grid(dsets, K = c(Inf, 30, 10, 5)) %>%
  mutate(out = pmap(list(N, R, E, seed, K), function(N, R, E, seed, K) {
    fit_sigma(simulate_oi(N = N, R = R, E = E, seed = seed), P = 2, K = K)
  })) %>%
  unnest(out)
saveRDS(res, "develop/20261124_circle_rank_sigma_map.rds")
