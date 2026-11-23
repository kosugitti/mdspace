# 円モデル専用の Stan(circle_rank.stan: 回転を足場で止める・使わない母数を置かない)を確かめる ---------------
# 20261122_ellipse_rank_std.R(閾値を潜在変数の広がりに合わせた版)では，3 次元のブリッジサンプリングの
# 繰り返しの SD が 0.25〜5.8 と大きかった。円モデルは回転しても値が変わらず，事後分布が弧の形になるため。
# circle_rank.stan は確率モデルを変えず，母数の作りだけを変えた(点 1 を原点・点 2〜N を下三角の行列)。
# 同じ 5 つのデータで，(1) 周辺尤度が前の版と同じ値になるか(体積の定数の確かめ)，
# (2) 3 次元のブリッジサンプリングの SD が小さくなるか，(3) 座標の有効標本数が増えるかを見る

pacman::p_load(tidyverse, cmdstanr, posterior, bridgesampling)
mod_circle <- cmdstanr::cmdstan_model("develop/circle_rank.stan", compile_model_methods = TRUE,
                                      force_recompile = TRUE)

# 原典の条件のスクリプトから simulate_oi() だけを読み込む ---------------------------------------------
local({
  ex <- parse("develop/20261122_ellipse_rank_recovery_oi1987.R")
  is_sim <- vapply(ex, function(e) {
    is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("simulate_oi"))
  }, logical(1))
  eval(ex[is_sim][[1]], envir = globalenv())
})

# 足場の P+1 個の点を選ぶ(対称部の Torgerson 解で，互いに離れ一つの超平面に乗らないもの) -----------------
pick_scaffold <- function(Z) {
  P <- ncol(Z)
  Zc <- sweep(Z, 2, colMeans(Z))
  idx <- which.max(rowSums(Zc^2))                  # 点 1: 重心から最も遠い
  R <- sweep(Z, 2, Z[idx, ])                       # 点 1 から見た位置
  for (p in seq_len(P)) {
    nxt <- which.max(rowSums(R^2))                 # これまでの点が張る空間から最も離れた点
    idx <- c(idx, nxt)
    u <- R[nxt, ] / sqrt(sum(R[nxt, ]^2))
    R <- R - tcrossprod(R %*% u, u)
  }
  idx
}

# 布置を足場の形(点 1 が原点・点 2〜N が下三角で対角が正)に回す -----------------------------------------
to_scaffold <- function(Z) {
  A <- sweep(Z[-1, , drop = FALSE], 2, Z[1, ])
  R <- qr.R(qr(t(A)))                              # t(A) = Q R なので A = t(R) t(Q)
  Lx <- t(R)
  sweep(Lx, 2, sign(diag(Lx[seq_len(ncol(Z)), , drop = FALSE])), `*`)
}

# データを Stan に渡す形にする(足場の点を先頭に並べ替える) ----------------------------------------------
prep_circle <- function(d, P) {
  d <- as.matrix(d)
  diag(d) <- NA
  N <- nrow(d)
  D0 <- d
  D0[is.na(D0)] <- t(D0)[is.na(D0)]
  diag(D0) <- 0
  Z0 <- stats::cmdscale((D0 + t(D0)) / 2, k = P)
  sc <- pick_scaffold(Z0)
  ord <- c(sc, setdiff(seq_len(N), sc))
  d <- d[ord, ord]
  Z0 <- Z0[ord, , drop = FALSE]
  use <- row(d) != col(d) & !is.na(d)
  vals <- sort(unique(d[use]))
  dat <- list(N = N, P = P, L = sum(use), row = row(d)[use], col = col(d)[use],
              K = length(vals), y = match(d[use], vals))
  make_init <- function(chain_id) {
    set.seed(chain_id)
    Z <- Z0 + matrix(rnorm(N * P, sd = 0.05), N, P)
    Z <- sweep(Z, 2, colMeans(Z))
    Z <- Z * sqrt(N / sum(Z^2))
    list(Lx = to_scaffold(Z), c_raw = rep(0, N), m_c = -2, s_c = 0.5, sig = 0.5,
         g = qnorm(seq_len(dat$K - 1) / dat$K))
  }
  list(dat = dat, make_init = make_init, ord = ord)
}

# 当てはめて周辺尤度を出す ---------------------------------------------------------------------------
fit_circle <- function(sim, P) {
  pr <- prep_circle(sim$d, P)
  D_true <- sim$D[pr$ord, pr$ord]
  od <- tempfile("fit_")
  dir.create(od)
  f <- mod_circle$sample(data = pr$dat, init = pr$make_init, chains = 4, parallel_chains = 4,
                         iter_warmup = 1000, iter_sampling = 2500, seed = 1, refresh = 0,
                         show_messages = FALSE, output_dir = od, sig_figs = 18)
  N <- pr$dat$N
  mu <- unclass(as.matrix(f$draws("mu", format = "matrix")))
  idx <- round(seq(1, nrow(mu), length.out = 1000))
  ktau <- vapply(idx, function(i) cor(mu[i, ], pr$dat$y, method = "kendall"), numeric(1))
  tau <- unclass(as.matrix(f$draws("tau", format = "matrix")))
  Dhat <- colMeans(t(apply(tau, 1, function(v) c(dist(matrix(v, N, P))))))
  ll <- rowSums(unclass(as.matrix(f$draws("log_lik", format = "matrix"))))
  sm <- f$summary("tau")
  rad <- exp(colMeans(unclass(as.matrix(f$draws("c", format = "matrix")))))
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
  tibble(logml = mean(b$logml),   # 回転と鏡映の体積は Stan の中で足してあるので補正は要らない
         logml_sd = sd(b$logml),
         loglik = mean(ll),
         divergent = sum(f$diagnostic_summary(quiet = TRUE)$num_divergent),
         tau_rhat_max = max(sm$rhat),
         tau_ess_min = min(sm$ess_bulk),
         sigma = mean(f$draws("sig")),
         ktau = mean(ktau),
         cor_dist = cor(Dhat, c(D_true[lower.tri(D_true)])),
         cor_r = cor(rad, sim$r[pr$ord]))   # 半径の回復(リカバリ用に後から足した。この下の 5 データの rds には無い)
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
    map_dfr(1:3, function(P) bind_cols(tibble(P = P), fit_circle(sim, P)))
  })) %>%
  unnest(out)
saveRDS(res, "develop/20261123_circle_rank_scaffold.rds")
