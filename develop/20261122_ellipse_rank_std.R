# 閾値を潜在変数の広がりで標準化した版(ellipse_rank_std.stan)で，次元の選び方を確かめる -----------------
# 20261122_ellipse_rank_cut.R では，閾値を mu の単位に置いても，事前分布の SD(1・2・4)で sigma が 0.5・1・3 と動き，
# SD 1 のときだけ 5 つすべてで真の 2 次元を選んだ。閾値の広がりは潜在変数 zeta = mu + 誤差 の広がりとほぼ同じになるので，
# 閾値の事前分布の幅を決めることが sigma を決めることになっていた。
# この版は閾値を cut = mean(mu) + sqrt(var(mu) + sig^2) * g，g は順序づけた N(0, 1) とし，人が決める幅が無い。
# 同じ 5 つのデータで，2 次元を選ぶか・sigma が真の誤差に近いかを見る

pacman::p_load(tidyverse, cmdstanr, posterior, bridgesampling)
source("develop/20261121_ellipse_rank_datasets.R")   # prep_data()
mod_std <- cmdstanr::cmdstan_model("develop/ellipse_rank_std.stan", compile_model_methods = TRUE,
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
fit_std <- function(sim, P) {
  pr <- prep_data(sim$d, "dissimilarity", P = P)
  # 初期値: 閾値 g は標準正規分布の分位点
  make_init <- function(chain_id) {
    ini <- pr$make_init(chain_id)
    ini$g <- qnorm(seq_len(pr$dat$K - 1) / pr$dat$K)
    ini$gam <- NULL
    ini
  }
  od <- tempfile("fit_")
  dir.create(od)
  f <- mod_std$sample(data = modifyList(pr$dat, list(use_size = 1, use_shape = 0)),
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
    map_dfr(1:3, function(P) bind_cols(tibble(P = P), fit_std(sim, P)))
  })) %>%
  unnest(out)
saveRDS(res, "develop/20261122_ellipse_rank_std.rds")
