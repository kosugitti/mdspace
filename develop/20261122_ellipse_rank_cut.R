# 閾値を mu の単位に置いた版(ellipse_rank_cut.stan)で，閾値の事前分布への敏感さを確かめ直す -------------
# 20261122_ellipse_rank_prior_loo.R では，閾値を sig の単位に置いた本体で，閾値の事前分布の SD(2・5・10)に
# よって sigma が 1.0・0.44・0.32 と動き，周辺尤度の選ぶ次元が入れ替わった。
# 原因は閾値の事前分布(0 が中心・sig の単位)が sig の強い事前分布として働いたこと。
# 閾値を mu の単位に置いて sig と切り離した版で，同じ 5 つのデータに当て，事前分布の SD(1・2・4)を変えても
# sigma と選ぶ次元が動かないかを見る。事前分布の平均は mu のおよその中心 1.3(布置は平方和 N に正規化している)

pacman::p_load(tidyverse, cmdstanr, posterior, bridgesampling)
source("develop/20261121_ellipse_rank_datasets.R")   # prep_data()
mod_cut <- cmdstanr::cmdstan_model("develop/ellipse_rank_cut.stan", compile_model_methods = TRUE,
                                   force_recompile = TRUE)

# 原典の条件のスクリプトから simulate_oi() だけを読み込む ---------------------------------------------
local({
  ex <- parse("develop/20261122_ellipse_rank_recovery_oi1987.R")
  is_sim <- vapply(ex, function(e) {
    is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("simulate_oi"))
  }, logical(1))
  eval(ex[is_sim][[1]], envir = globalenv())
})

# 当てはめて周辺尤度を出す ---------------------------------------------------------------------------
fit_cut <- function(sim, P, cut_sd, cut_m = 1.3) {
  pr <- prep_data(sim$d, "dissimilarity", P = P)
  # 初期値: 本体の閾値 gam(sig の単位)を mu の単位に直す
  make_init <- function(chain_id) {
    ini <- pr$make_init(chain_id)
    ini$cut <- ini$gam * ini$sig
    ini$gam <- NULL
    ini
  }
  od <- tempfile("fit_")
  dir.create(od)
  f <- mod_cut$sample(data = modifyList(pr$dat, list(use_size = 1, use_shape = 0,
                                                     cut_m = cut_m, cut_sd = cut_sd)),
                      init = make_init, chains = 4, parallel_chains = 4, iter_warmup = 1000,
                      iter_sampling = 2500, seed = 1, refresh = 0, show_messages = FALSE,
                      output_dir = od, sig_figs = 18)
  N <- pr$dat$N
  mu <- unclass(as.matrix(f$draws("mu", format = "matrix")))
  idx <- round(seq(1, nrow(mu), length.out = 1000))
  ktau <- vapply(idx, function(i) cor(mu[i, ], pr$dat$y, method = "kendall"), numeric(1))
  tau <- unclass(as.matrix(f$draws("tau", format = "matrix")))
  Dhat <- colMeans(t(apply(tau, 1, function(v) c(dist(matrix(v, N, P))))))
  ll <- rowSums(unclass(as.matrix(f$draws("log_lik", format = "matrix"))))
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
         loglik = mean(ll),
         divergent = sum(f$diagnostic_summary(quiet = TRUE)$num_divergent),
         sigma = mean(f$draws("sig")),
         ktau = mean(ktau),
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
    expand_grid(cut_sd = c(1, 2, 4), P = 1:3) %>%
      mutate(fit = map2(P, cut_sd, function(P, s) fit_cut(sim, P, s))) %>%
      unnest(fit)
  })) %>%
  unnest(out)
saveRDS(res, "develop/20261122_ellipse_rank_cut.rds")
