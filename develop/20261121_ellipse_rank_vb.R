# 楕円モデルの順位尤度版を，最適化(事後確率最大値)・Pathfinder・変分推論(ADVI)で推定して MCMC と突き合わせる -----
# データの準備と初期値は develop/20261121_ellipse_rank_try.R と同じ

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

# 推定 -----------------------------------------------------------------------------
mod <- cmdstanr::cmdstan_model("develop/ellipse_rank.stan")
mcmc <- readRDS("develop/20261121_ellipse_rank_fit.rds")

# 最適化: 初期値を30本変え，対数事後密度(ヤコビアンなし)の最大を取る
opts <- map(1:30, function(i) {
  tryCatch(mod$optimize(data = dat, init = list(make_init(i)), seed = i, jacobian = FALSE,
                        iter = 10000, refresh = 0, show_messages = FALSE),
           error = function(e) NULL)
}) %>% compact()
lps <- map_dbl(opts, ~ .x$lp())
opt <- opts[[which.max(lps)]]

pf <- mod$pathfinder(data = dat, init = map(1:4, make_init), num_paths = 4, seed = 1,
                     refresh = 0, show_messages = FALSE)
vb <- mod$variational(data = dat, init = list(make_init(1)), seed = 1, iter = 20000,
                      refresh = 0, show_messages = FALSE)

# 比べる量: 超母数と，点間距離(向きによらない) ------------------------------------------------
dist_draws <- function(tau_mat) {
  # 各行が一つのサンプルの tau(N x P を列優先で並べたもの)
  t(apply(tau_mat, 1, function(v) c(dist(matrix(v, N, P)))))
}
pars <- c("m_c", "s_c", "s_e", "sig")
tab <- bind_rows(
  tibble(method = "MCMC", par = pars, value = map_dbl(pars, ~ mean(mcmc$draws(.x)))),
  tibble(method = "Pathfinder", par = pars, value = map_dbl(pars, ~ mean(pf$draws(.x)))),
  tibble(method = "ADVI", par = pars, value = map_dbl(pars, ~ mean(vb$draws(.x)))),
  tibble(method = "Optimize", par = pars, value = map_dbl(pars, ~ opt$mle(.x)))
) %>% pivot_wider(names_from = method, values_from = value)

D_mcmc <- colMeans(dist_draws(mcmc$draws("tau", format = "matrix")))
D_pf <- colMeans(dist_draws(pf$draws("tau", format = "matrix")))
D_vb <- colMeans(dist_draws(vb$draws("tau", format = "matrix")))
D_opt <- c(dist(matrix(opt$mle("tau"), N, P)))
cor_D <- cor(cbind(MCMC = D_mcmc, Pathfinder = D_pf, ADVI = D_vb, Optimize = D_opt))

saveRDS(list(tab = tab, cor_D = cor_D, lps = lps, opt_phi = matrix(opt$mle("phi"), N, P),
             opt_tau = matrix(opt$mle("tau"), N, P), lab = lab),
        "develop/20261121_ellipse_rank_vb.rds")
