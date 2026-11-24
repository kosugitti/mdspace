# 円モデルの順位尤度版(circle_rank.stan)の事後確率最大値: 布置の初期値を乱数にして局所解の数を見る -----------
# 20261123_circle_rank_approx.R では，Torgerson 解の近くから始めた最適化は 2〜3 次元でほぼ同じ山に落ちた。
# Stan の乱数の初期値は，131 個の閾値の間隔まで乱数になって確率 0 のセルが出るので 1 歩目から計算できなかった。
# ここでは閾値・sig・半径の初期値はそのままにして，布置だけを乱数(独立な標準正規・重心 0・平方和 N)にする。
# R で書くときに「Torgerson 解から 1 本」で足りるかどうかを決める材料
# 1 次元は点が重なる所の角で勾配法が止まるので(20261123_circle_rank_approx.R の注)，2・3 次元だけ

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

# 1 つのデータ・次元で，Torgerson 解から 1 本と乱数の布置から n_random 本 --------------------------------
run_starts <- function(sim, P, n_random = 50) {
  pr <- prep_circle(sim$d, P)
  N <- pr$dat$N
  K <- pr$dat$K
  od <- tempfile("fit_")
  dir.create(od)
  D_true <- sim$D[pr$ord, pr$ord][lower.tri(diag(N))]
  r_true <- sim$r[pr$ord]
  make_init <- function(i) {
    if (i == 0) {
      Lx <- pr$make_init(1)$Lx                       # Torgerson 解の近く
    } else {
      set.seed(1000 + i)
      Z <- matrix(rnorm(N * P), N, P)                # 乱数の布置
      Z <- sweep(Z, 2, colMeans(Z))
      Lx <- to_scaffold(Z * sqrt(N / sum(Z^2)))
    }
    list(Lx = Lx, a_raw = rep(0, N), s_a = 0.3, sig = 0.5, g = qnorm(seq_len(K - 1) / K))
  }
  map_dfr(0:n_random, function(i) {
    o <- tryCatch(mod$optimize(data = pr$dat, init = list(make_init(i)), seed = 1, jacobian = FALSE,
                               iter = 20000, refresh = 0, show_messages = FALSE, output_dir = od),
                  error = function(e) NULL)
    lp <- if (is.null(o)) NA_real_ else tryCatch(o$lp(), error = function(e) NA_real_)
    if (is.na(lp)) {
      return(tibble(start = i, lp = NA_real_, D = list(NULL), sigma = NA_real_,
                    cor_D_true = NA_real_, cor_r_true = NA_real_))
    }
    Dh <- c(dist(matrix(o$mle("tau"), N, P)))
    tibble(start = i, lp = lp, D = list(Dh), sigma = unname(o$mle("sig")),
           cor_D_true = cor(Dh, D_true), cor_r_true = cor(o$mle("r"), r_true))
  })
}

# データと条件(20261123_circle_rank_scaffold.R と同じ 5 つ) ----------------------------------------
dsets <- tribble(
  ~label,               ~N, ~R,  ~E,  ~seed,
  "N12・E0.5・データ1", 12, 0.3, 0.5, 2001,
  "N12・E0.5・データ2", 12, 0.3, 0.5, 2002,
  "N10・E0.7・データ1", 10, 0.3, 0.7, 2009,
  "N10・E0.7・データ2", 10, 0.3, 0.7, 2010,
  "N12・E0.3(対照)",   12, 0.3, 0.3, 1007
)
res <- expand_grid(dsets, P = 2:3) %>%
  mutate(out = pmap(list(N, R, E, seed, P), function(N, R, E, seed, P) {
    run_starts(simulate_oi(N = N, R = R, E = E, seed = seed), P)
  })) %>%
  unnest(out) %>%
  group_by(label, P) %>%
  mutate(dlp = max(lp, na.rm = TRUE) - lp,           # 最良の解との対数事後密度の差
         cor_best = map_dbl(D, function(d) {         # 最良の解との点間距離の相関
           if (is.null(d)) NA_real_ else cor(d, D[[which.max(lp)]])
         })) %>%
  ungroup() %>%
  select(-D)
saveRDS(res, "develop/20261124_circle_rank_map_starts.rds")
