# 次元の選び方: 閾値の事前分布への敏感さと LOO ------------------------------------------------------
# 20261122_ellipse_rank_recovery_harsh.R では，E = 0.5〜1.0 で周辺尤度が真の 2 次元より 1 次元を選んだ
# (点間距離の相関は 2 次元のほうが高いのに)。罰則の出どころを探る
#   (1) 閾値の事前分布の SD(2, 5, 10)で周辺尤度の次元の差が動くか
#       座標の事前分布の SD は確かめない: 布置は平方和 N に正規化するので，等方の正規分布の SD は
#       正規化後の布置の分布を変えない(基準の対象の制約も尺度に依らない)
#   (2) LOO(セルごとの予測)で比べると何次元を選ぶか
# データ: 厳しい条件で 1 次元を選んだ N = 12・E = 0.5 の 2 つと N = 10・E = 0.7 の 2 つ，
#         対照に原典の条件(N = 12・R = 0.3・E = 0.3)で 2 次元を選んだ 1 つ

pacman::p_load(tidyverse, cmdstanr, posterior, bridgesampling, loo)
source("develop/20261121_ellipse_rank_datasets.R")   # prep_data()
mod_gsd <- cmdstanr::cmdstan_model("develop/ellipse_rank_gsd.stan", compile_model_methods = TRUE,
                                   force_recompile = TRUE)

# 原典の条件のスクリプトから simulate_oi() だけを読み込む ---------------------------------------------
local({
  ex <- parse("develop/20261122_ellipse_rank_recovery_oi1987.R")
  is_sim <- vapply(ex, function(e) {
    is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("simulate_oi"))
  }, logical(1))
  eval(ex[is_sim][[1]], envir = globalenv())
})

# 当てはめて周辺尤度と LOO を出す ---------------------------------------------------------------------
fit_gsd <- function(sim, P, gam_sd) {
  pr <- prep_data(sim$d, "dissimilarity", P = P)
  od <- tempfile("fit_")
  dir.create(od)
  f <- mod_gsd$sample(data = modifyList(pr$dat, list(use_size = 1, use_shape = 0, gam_sd = gam_sd)),
                      init = pr$make_init, chains = 4, parallel_chains = 4, iter_warmup = 1000,
                      iter_sampling = 2500, seed = 1, refresh = 0, show_messages = FALSE,
                      output_dir = od, sig_figs = 18)
  N <- pr$dat$N
  tau <- unclass(as.matrix(f$draws("tau", format = "matrix")))
  Dhat <- colMeans(t(apply(tau, 1, function(v) c(dist(matrix(v, N, P))))))
  ll <- f$draws("log_lik", format = "draws_array")
  lo <- loo::loo(ll, r_eff = loo::relative_eff(exp(ll)))
  f$init_model_methods()
  U <- unclass(as.matrix(f$unconstrain_draws(format = "draws_matrix")))
  colnames(U) <- paste0("u", seq_len(ncol(U)))
  lb <- setNames(rep(-Inf, ncol(U)), colnames(U))
  ub <- setNames(rep(Inf, ncol(U)), colnames(U))
  # 提案の点で sig が 0 に近いと mu / sig が inf になり例外が出るので，そのときは -Inf を返す
  lp <- function(x, data) {
    tryCatch(f$log_prob(unconstrained_variables = x, jacobian = TRUE), error = function(e) -Inf)
  }
  b <- bridge_sampler(samples = U, log_posterior = lp, data = NULL, lb = lb, ub = ub,
                      method = "warp3", repetitions = 3, silent = TRUE)
  tibble(logml = mean(b$logml) + P * log(2),   # 鏡映(基準の対象を正にした半空間)の補正
         logml_sd = sd(b$logml),
         elpd_loo = lo$estimates["elpd_loo", "Estimate"],
         se_elpd_loo = lo$estimates["elpd_loo", "SE"],
         p_loo = lo$estimates["p_loo", "Estimate"],
         n_khat_bad = sum(lo$diagnostics$pareto_k > 0.7),
         elpd_pointwise = list(lo$pointwise[, "elpd_loo"]),
         divergent = sum(f$diagnostic_summary(quiet = TRUE)$num_divergent),
         sigma = mean(f$draws("sig")),
         cor_dist = cor(Dhat, c(sim$D[lower.tri(sim$D)])))
}

# データと条件 --------------------------------------------------------------------------------------
dsets <- tribble(
  ~label,               ~N, ~R,  ~E,  ~seed,
  "N12・E0.5・データ1", 12, 0.3, 0.5, 2001,
  "N12・E0.5・データ2", 12, 0.3, 0.5, 2002,
  "N10・E0.7・データ1", 10, 0.3, 0.7, 2009,
  "N10・E0.7・データ2", 10, 0.3, 0.7, 2010,
  "N12・E0.3(対照)",   12, 0.3, 0.3, 1007
)
res <- dsets %>%
  mutate(out = pmap(list(N, R, E, seed), function(N, R, E, seed) {
    sim <- simulate_oi(N = N, R = R, E = E, seed = seed)
    expand_grid(gam_sd = c(2, 5, 10), P = 1:3) %>%
      mutate(fit = map2(P, gam_sd, function(P, g) fit_gsd(sim, P, g))) %>%
      unnest(fit)
  })) %>%
  unnest(out)
saveRDS(res, "develop/20261122_ellipse_rank_prior_loo.rds")
