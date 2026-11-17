# circle_mds() のパラメタリカバリ(真の布置と円から作ったデータを当て直す) ------------------
# mdspace のプロジェクトのフォルダを作業ディレクトリにして実行する
# 判定:
#   局所解 = 推定の S が「真の布置で計算した S」より大きい(もっと良い解があるのに着けていない)
#   回復   = 正規化と直交回転(鏡映を含む)をそろえた座標の誤差，半径の相関・誤差
pacman::p_load(tidyverse, pkgload, parallel)

pkgload::load_all(".")

# 条件 ------------------------------------------------------------------------
# 次元は対象の数より小さくないと解けないので，n = 5 は 2・3次元，n = 10 は 2・3・5次元
conds <- tidyr::expand_grid(n = c(5, 10, 20), p = c(2, 3, 5, 10),
                            noise = c(0, 0.1)) %>%
  dplyr::filter(p < n)
nrep <- 20

# 真の布置と円からデータを作る ---------------------------------------------------------
# 座標は標準正規，半径は 0 から距離の中央値の 0.3 倍までの一様乱数
# 類似度 = -(m_jk + 誤差)。誤差の SD は m の SD の noise 倍
make_data <- function(n, p, noise) {
  X <- matrix(stats::rnorm(n * p), n, p)
  D <- as.matrix(stats::dist(X))
  r <- stats::runif(n, 0, 0.3 * stats::median(D[upper.tri(D)]))
  M <- D - outer(r, r, "-")
  off <- row(M) != col(M)
  E <- matrix(stats::rnorm(n * n, 0, noise * stats::sd(M[off])), n, n)
  S <- -(M + E)
  diag(S) <- NA
  lab <- paste0("o", seq_len(n))
  dimnames(S) <- list(lab, lab)
  return(list(S = S, X = X, r = r))
}

# 真の布置で S を計算する(推定と同じ正規化・同じ単調回帰) ------------------------------------
true_stress <- function(dat) {
  S <- dat$S
  use <- row(S) != col(S)
  cells <- list(s = S[use], j = row(S)[use], k = col(S)[use], n = nrow(S),
                decreasing = TRUE)
  return(mdspace:::.circle_stress(c(dat$X, dat$r), cells, ncol(dat$X)))
}

# 回復の指標 ------------------------------------------------------------------------
recovery <- function(fit, dat) {
  tr <- mdspace:::.circle_normalize(dat$X, dat$r)
  sv <- svd(t(fit$conf) %*% tr$conf)
  rot <- fit$conf %*% sv$u %*% t(sv$v)
  return(tibble::tibble(
    conf_rmse = sqrt(mean((rot - tr$conf)^2)),
    conf_max = max(abs(rot - tr$conf)),
    radius_cor = suppressWarnings(stats::cor(fit$radius, tr$radius)),
    radius_max = max(abs(fit$radius - tr$radius))
  ))
}

# 1回分 -------------------------------------------------------------------------
one_run <- function(i, n, p, noise) {
  set.seed(1000 * n + 100 * p + 10 * noise * 10 + i)
  dat <- make_data(n, p, noise)
  s_true <- true_stress(dat)
  t0 <- proc.time()[["elapsed"]]
  fit <- circle_mds(dat$S, ndim = p, nrandom = 10, seed = i)
  sec <- proc.time()[["elapsed"]] - t0
  return(dplyr::bind_cols(
    tibble::tibble(n = n, p = p, noise = noise, rep = i, s_true = s_true,
                   s_fit = fit$stress, s_rational = fit$starts[1],
                   n_starts_at_best = sum(fit$starts <= fit$stress + 1e-4),
                   best_start = fit$best_start, sec = sec),
    recovery(fit, dat)
  ))
}

# 実行(条件 x 反復を並列) ------------------------------------------------------------
jobs <- conds %>% tidyr::crossing(rep = seq_len(nrep))
res <- parallel::mclapply(seq_len(nrow(jobs)), function(k) {
  one_run(jobs$rep[k], jobs$n[k], jobs$p[k], jobs$noise[k])
}, mc.cores = max(1, parallel::detectCores() - 2)) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(local_min = s_fit > s_true + 1e-4,
                rational_local = s_rational > s_true + 1e-4)

dir.create("develop/sim_results", showWarnings = FALSE)
readr::write_csv(res, "develop/sim_results/20261117_circle_recovery.csv")

# 条件ごとのまとめ ----------------------------------------------------------------------
summary_tab <- res %>%
  dplyr::group_by(n, p, noise) %>%
  dplyr::summarise(
    local_min = sum(local_min),
    rational_local = sum(rational_local),
    s_true = mean(s_true),
    s_fit = mean(s_fit),
    conf_rmse = median(conf_rmse),
    conf_max = max(conf_max),
    radius_cor = median(radius_cor),
    sec = mean(sec),
    .groups = "drop"
  )
readr::write_csv(summary_tab, "develop/sim_results/20261117_circle_recovery_summary.csv")
summary_tab
