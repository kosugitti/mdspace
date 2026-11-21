# 楕円モデルの順位尤度版(develop/ellipse_rank.stan・正規の尺度＋母数の空間の制約)を職業データ(書籍表1.19)で試す ----
# データは類似度なので符号を反転して非類似度の順位にする(順序尺度として扱うので値は使わない)

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

# MCMC ------------------------------------------------------------------------------
mod <- cmdstanr::cmdstan_model("develop/ellipse_rank.stan")
fit <- mod$sample(data = dat, init = make_init, chains = 4, parallel_chains = 4,
                  iter_warmup = 1000, iter_sampling = 1000, seed = 1, refresh = 0,
                  show_messages = FALSE)

# 診断 ------------------------------------------------------------------------------
diag_sum <- fit$diagnostic_summary()
summ <- fit$summary(c("m_c", "s_c", "s_e", "sig", "tau", "phi", "c"))
saveRDS(list(fit_summary = summ, diag = diag_sum, dat = dat, lab = lab),
        "develop/20261121_ellipse_rank_try.rds")
fit$save_object("develop/20261121_ellipse_rank_fit.rds")
