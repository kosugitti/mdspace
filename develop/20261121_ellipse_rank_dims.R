# 次元数の比較と適合度: 4つのデータ × P = 1, 2, 3, 4 × 円・楕円 ----------------------------------------
# 周辺尤度(ブリッジサンプリング)で次元とモデルを比べ，適合度として sigma(布置の広がりに比べた誤差)と
# 順位の一致率(mu とデータの Kendall の tau-b の事後分布)を出す。
# 軸の入れ替えの対称性で事後分布に最大 P! 個の山があり，鎖が一つの山だけを回ると周辺尤度は最大 log(P!)
# 過小に見積もられる(基準の対象を軸ごとに選んでいるので山は厳密には同じ形ではない)。補正前と補正後の両方を出す
# (2026-11-22 訂正) 鏡映: 基準の対象で軸ごとに半分に切っているので P log 2 を足す。軸の入れ替え log(P!) は楕円だけ。
# 円モデルは回転が自由で鎖が回転の方向を回るので入れ替えの補正は要らない
# P = 1 では楕円と円は同じモデルなので円だけを当てる

pacman::p_load(tidyverse, cmdstanr, posterior, bridgesampling)
source("develop/20261121_ellipse_rank_datasets.R")   # prep_data()・harshman()・mod

datasets <- list()
datasets[["職業"]] <- as.matrix(read.csv("tests/testthat/fixtures/okada_imaizumi_2024_table1_19.csv",
                                         row.names = 1, check.names = FALSE))
x <- read.csv("tests/testthat/fixtures/okada_imaizumi_2024_table1_14.csv", check.names = FALSE)
for (b in split(x, x$source)) {
  m <- as.matrix(b[, -(1:2)])
  rownames(m) <- b$from
  datasets[[paste0("国の印象", b$source[1])]] <- m
}
data(carswitching, package = "asymmetry")
datasets[["車の乗り換え"]] <- harshman(as.matrix(carswitching))$rescaleddata

fit_dim <- function(s, name, P, use_shape) {
  pr <- prep_data(s, "similarity", P = P)
  od <- tempfile("fit_")
  dir.create(od)
  f <- mod$sample(data = modifyList(pr$dat, list(use_size = 1, use_shape = use_shape)),
                  init = pr$make_init, chains = 4, parallel_chains = 4, iter_warmup = 1000,
                  iter_sampling = 2500, seed = 1, refresh = 0, show_messages = FALSE,
                  output_dir = od, sig_figs = 18)   # 和が 0 のベクトルを CSV の丸めで壊さない(P = 3 で必要)
  # 順位の一致率: サンプルを間引いて(1,000本)Kendall の tau-b を計算する(ブリッジサンプリングの前に)
  mu <- unclass(as.matrix(f$draws("mu", format = "matrix")))
  idx <- round(seq(1, nrow(mu), length.out = 1000))
  ktau <- vapply(idx, function(i) cor(mu[i, ], pr$dat$y, method = "kendall"), numeric(1))
  f$init_model_methods()
  U <- unclass(as.matrix(f$unconstrain_draws(format = "draws_matrix")))
  colnames(U) <- paste0("u", seq_len(ncol(U)))
  lb <- setNames(rep(-Inf, ncol(U)), colnames(U))
  ub <- setNames(rep(Inf, ncol(U)), colnames(U))
  lp <- function(x, data) f$log_prob(unconstrained_variables = x, jacobian = TRUE)
  b <- bridge_sampler(samples = U, log_posterior = lp, data = NULL, lb = lb, ub = ub,
                      method = "warp3", repetitions = 5, silent = TRUE)
  sig <- as.vector(f$draws("sig"))
  tibble(data = name, P = P, model = if (use_shape == 1) "楕円" else "円",
         logml = mean(b$logml), logml_sd = sd(b$logml), logml_corr = mean(b$logml) + P * log(2) + if (use_shape == 1) lfactorial(P) else 0,
         divergent = sum(f$diagnostic_summary(quiet = TRUE)$num_divergent),
         sigma = mean(sig), sigma_lo = quantile(sig, 0.025), sigma_hi = quantile(sig, 0.975),
         ktau_lo = quantile(ktau, 0.025), ktau_hi = quantile(ktau, 0.975), ktau = mean(ktau))
}

grid <- expand_grid(name = names(datasets), P = 1:4, use_shape = 0:1) %>%
  filter(!(P == 1 & use_shape == 1))
res <- pmap_dfr(grid, function(name, P, use_shape) fit_dim(datasets[[name]], name, P, use_shape))
saveRDS(res, "develop/20261121_ellipse_rank_dims.rds")
# (注) 11/21 の実行では ktau_lo・ktau_hi が平均と同じ値になった。tibble() は引数を順に評価し，先に作った列 ktau(平均)を
# 後の quantile() が参照していたため(平均値そのものは正しい)。区間を先に作るよう並べ替えた
