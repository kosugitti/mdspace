# 下位モデル(楕円・円・等面積の楕円・対称)の周辺尤度をブリッジサンプリングで求めて比べる ---------------
# LOO は閾値がセルごとにほぼ1個ずつあるため崩れる(p_loo がセル数を超える)ので，周辺尤度で比べる
# 4モデルは母数と事前分布の並びが同じ(使わない散らばりも事前分布だけで動く)ので，~ で落ちる定数は比で打ち消し合う

pacman::p_load(tidyverse, cmdstanr, posterior)

# データ -----------------------------------------------------------------------------
s <- as.matrix(read.csv("tests/testthat/fixtures/okada_imaizumi_2024_table1_19.csv",
                        row.names = 1, check.names = FALSE))
lab <- rownames(s)
N <- nrow(s)
P <- 2
use <- row(s) != col(s) & !is.na(s)
d <- -s[use]                                   # 非類似度(順位だけを使う)
vals <- sort(unique(d))
y <- match(d, vals)                            # 下から何番目の値か(タイは同じ番号)

# 基準の対象: 対称部の Torgerson 解で，軸ごとに重心から最も遠い対象 ------------------------------
D0 <- max(s, na.rm = TRUE) - s
diag(D0) <- 0
Z0 <- stats::cmdscale((D0 + t(D0)) / 2, k = P)
anchor <- apply(abs(sweep(Z0, 2, colMeans(Z0))), 2, which.max)
# 基準の対象の座標が正になるよう軸の向きをそろえる(初期値用)
Z0 <- sweep(Z0, 2, sign(Z0[cbind(anchor, seq_len(P))]), `*`)

dat <- list(N = N, P = P, L = sum(use), row = row(s)[use], col = col(s)[use],
            K = length(vals), y = y, anchor = anchor, use_size = 1, use_shape = 1)

# 初期値: Torgerson 解・小さな円(mu が正になるように)・閾値は初期の log(mu) に沿って並べる -----------------
make_init <- function(chain_id) {
  set.seed(chain_id)
  Z <- Z0 + matrix(rnorm(N * P, sd = 0.05), N, P)
  Z <- sweep(Z, 2, colMeans(Z))
  Z <- Z * sqrt(N / sum(Z^2))
  yr <- matrix(0, N - 1, P)
  v <- numeric(P)
  for (p in seq_len(P)) {
    others <- setdiff(seq_len(N), anchor[p])
    yr[, p] <- Z[others, p]
    v[p] <- log(Z[anchor[p], p] - mean(Z[others, p]))
  }
  dl <- sort(sqrt(rowSums((Z[dat$row, ] - Z[dat$col, ])^2)))
  sig <- 0.5
  gam <- seq(dl[2], dl[length(dl) - 1], length.out = dat$K - 1) / sig
  list(yr = yr, v = v, c_raw = rep(0, N), et = replicate(N, rep(0, P), simplify = FALSE),
       m_c = -2, s_c = 0.5, s_e = 0.5, sig = sig, gam = gam)
}

pacman::p_load(bridgesampling)
mod <- cmdstanr::cmdstan_model("develop/ellipse_rank.stan", compile_model_methods = TRUE,
                               force_recompile = TRUE)
specs <- tribble(
  ~model,         ~use_size, ~use_shape,
  "楕円(II)",          1,          1,
  "円(III)",           1,          0,
  "等面積の楕円",      0,          1,
  "対称",              0,          0
)
bs_one <- function(model, use_size, use_shape) {
  od <- tempfile("fit_")
  dir.create(od)
  f <- mod$sample(data = modifyList(dat, list(use_size = use_size, use_shape = use_shape)),
                  init = make_init, chains = 4, parallel_chains = 4, iter_warmup = 1000,
                  iter_sampling = 2500, seed = 1, refresh = 0, show_messages = FALSE,
                  output_dir = od)
  f$init_model_methods()
  # 制約のない空間のサンプル(行がサンプル)
  ud <- f$unconstrain_draws(format = "draws_matrix")
  U <- unclass(as.matrix(ud))
  colnames(U) <- paste0("u", seq_len(ncol(U)))
  lb <- setNames(rep(-Inf, ncol(U)), colnames(U))
  ub <- setNames(rep(Inf, ncol(U)), colnames(U))
  lp <- function(x, data) f$log_prob(unconstrained_variables = x, jacobian = TRUE)
  b <- bridge_sampler(samples = U, log_posterior = lp, data = NULL, lb = lb, ub = ub,
                      method = "warp3", repetitions = 5, silent = TRUE)
  div <- sum(f$diagnostic_summary(quiet = TRUE)$num_divergent)
  list(model = model, logml = b$logml, logml_sd = sd(b$logml), div = div)
}
res <- pmap(specs, bs_one)
tab <- map_dfr(res, ~ tibble(model = .x$model, logml = mean(.x$logml), logml_sd = .x$logml_sd,
                             divergent = .x$div)) %>%
  mutate(diff_from_best = logml - max(logml))
saveRDS(tab, "develop/20261121_ellipse_rank_bridge.rds")
