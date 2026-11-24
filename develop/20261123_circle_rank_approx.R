# 円モデルの順位尤度版(circle_rank.stan: 足場・標準化した閾値・半径は差だけ)を MCMC 以外で推定して突き合わせる ---
# 本命は R による近似計算(事後確率最大値)。その前に，Stan の最適化・Pathfinder・ADVI が MCMC と合うかを見る。
# 同じ 5 つのデータ・1〜3 次元で，(1) 点間距離・半径・sigma・s_a が MCMC と合うか，
# (2) 最適化の初期値(Torgerson 解の近く 10 本・Stan の乱数 10 本)が同じ山に落ちるかを確かめる。
# MCMC は adapt_delta = 0.99(閾値の下の裾の発散を避ける。推定値は 0.8 と変わらないことを確かめ済み)

pacman::p_load(tidyverse, cmdstanr, posterior)
mod <- cmdstanr::cmdstan_model("develop/circle_rank.stan")

# simulate_oi() と，足場のスクリプトの pick_scaffold()・to_scaffold()・prep_circle() を読み込む -------------
load_defs <- function(file, names) {
  for (e in parse(file)) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]]) %in% names) {
      eval(e, envir = globalenv())
    }
  }
}
load_defs("develop/20261122_ellipse_rank_recovery_oi1987.R", "simulate_oi")
load_defs("develop/20261123_circle_rank_scaffold.R", c("pick_scaffold", "to_scaffold", "prep_circle"))

# 推定値から比べる量を取り出す(draws は行が 1 つの標本) -------------------------------------------
summ <- function(tau, r, sig, s_a, N, P) {
  D <- t(apply(tau, 1, function(v) c(dist(matrix(v, N, P)))))
  list(D = colMeans(D), r = colMeans(r), sig = mean(sig), s_a = mean(s_a))
}
get_draws <- function(f, N, P) {
  m <- function(v) unclass(as.matrix(f$draws(v, format = "matrix")))
  summ(m("tau"), m("r"), m("sig"), m("s_a"), N, P)
}

# 1 つのデータ・次元で全部の推定法を回す ------------------------------------------------------------
fit_all <- function(sim, P) {
  pr <- prep_circle(sim$d, P)
  N <- pr$dat$N
  init <- function(chain_id) {
    i0 <- pr$make_init(chain_id)
    i0[c("c_raw", "m_c", "s_c")] <- NULL
    i0$a_raw <- rep(0, N)
    i0$s_a <- 0.3
    i0
  }
  od <- tempfile("fit_")
  dir.create(od)
  D_true <- sim$D[pr$ord, pr$ord][lower.tri(diag(N))]
  r_true <- sim$r[pr$ord]
  out <- list()

  # MCMC(基準)
  t0 <- Sys.time()
  f <- mod$sample(data = pr$dat, init = init, chains = 4, parallel_chains = 4, iter_warmup = 1000,
                  iter_sampling = 2500, seed = 1, refresh = 0, show_messages = FALSE,
                  adapt_delta = 0.99, output_dir = od)
  out$MCMC <- c(get_draws(f, N, P), sec = as.numeric(difftime(Sys.time(), t0, units = "secs")))

  # 最適化(事後確率最大値・ヤコビアンなし): Torgerson 解の近く 10 本と Stan の乱数の初期値 10 本
  t0 <- Sys.time()
  opts <- map(1:20, function(i) {
    ini <- if (i <= 10) list(init(i)) else NULL
    tryCatch(mod$optimize(data = pr$dat, init = ini, seed = i, jacobian = FALSE, iter = 20000,
                          refresh = 0, show_messages = FALSE, output_dir = od),
             error = function(e) NULL)
  })
  # 失敗した回は cmdstanr が例外を出さずに返し，lp() を呼んだときに初めて止まるので，ここで NA にする
  lps <- vapply(opts, function(o) {
    if (is.null(o)) return(NA_real_)
    tryCatch(o$lp(), error = function(e) NA_real_)
  }, numeric(1))
  ok <- !is.na(lps)
  # 1 次元では，点が重なる所(|x_j - x_k| の角)で直線探索が止まり，20 本とも失敗することがある(下の注)
  D_best <- NULL
  if (any(ok)) {
    best <- opts[[which.max(lps)]]
    mle <- function(v) matrix(best$mle(v), nrow = 1)
    out$MAP <- c(summ(mle("tau"), mle("r"), best$mle("sig"), best$mle("s_a"), N, P),
                 sec = as.numeric(difftime(Sys.time(), t0, units = "secs")) / 20)
    D_best <- out$MAP$D
  }
  # 最良の山からのずれ: 対数事後密度の差と，最良の解との点間距離の相関
  opt_tab <- tibble(start = rep(c("torgerson", "random"), each = 10), lp = lps,
                    cor_best = map_dbl(seq_len(20), function(i) {
                      if (!ok[i] || is.null(D_best)) return(NA_real_)
                      cor(c(dist(matrix(opts[[i]]$mle("tau"), N, P))), D_best)
                    }))

  # Pathfinder(4 経路)
  t0 <- Sys.time()
  pf <- tryCatch(mod$pathfinder(data = pr$dat, init = map(1:4, init), num_paths = 4, seed = 1,
                                refresh = 0, show_messages = FALSE, output_dir = od),
                 error = function(e) NULL)
  if (!is.null(pf)) {
    out$Pathfinder <- c(get_draws(pf, N, P), sec = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }

  # ADVI(平均場)
  t0 <- Sys.time()
  vb <- tryCatch(mod$variational(data = pr$dat, init = list(init(1)), seed = 1, iter = 20000,
                                 refresh = 0, show_messages = FALSE, output_dir = od),
                 error = function(e) NULL)
  if (!is.null(vb)) {
    out$ADVI <- c(get_draws(vb, N, P), sec = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }

  ref <- out$MCMC
  tab <- imap_dfr(out, function(o, nm) {
    tibble(method = nm,
           cor_D_mcmc = cor(o$D, ref$D),
           cor_r_mcmc = cor(o$r, ref$r),
           sigma = o$sig, s_a = o$s_a,
           r_max = max(o$r),
           cor_D_true = cor(o$D, D_true),
           cor_r_true = cor(o$r, r_true),
           sec = o$sec)
  })
  list(tab = tab, opt = opt_tab)
}

# 注(1 次元の最適化): N12・E0.3 の 1 次元では，L-BFGS が勾配のノルム約 19 のまま直線探索で止まった。
# 止まる直前の解では 3 つの点の間隔が 0.002 と 0.0007 で，点が一つに寄っていた。1 次元の距離 |x_j - x_k| は
# 点が重なる所で微分できないので，勾配法はそこで進めなくなる(1 次元の非計量 MDS で知られた性質)

# データと条件(20261123_circle_rank_scaffold.R と同じ 5 つ) ----------------------------------------
dsets <- tribble(
  ~label,               ~N, ~R,  ~E,  ~seed,
  "N12・E0.5・データ1", 12, 0.3, 0.5, 2001,
  "N12・E0.5・データ2", 12, 0.3, 0.5, 2002,
  "N10・E0.7・データ1", 10, 0.3, 0.7, 2009,
  "N10・E0.7・データ2", 10, 0.3, 0.7, 2010,
  "N12・E0.3(対照)",   12, 0.3, 0.3, 1007
)
res <- expand_grid(dsets, P = 1:3) %>%
  mutate(out = pmap(list(N, R, E, seed, P), function(N, R, E, seed, P) {
    fit_all(simulate_oi(N = N, R = R, E = E, seed = seed), P)
  }))
saveRDS(res, "develop/20261123_circle_rank_approx.rds")
