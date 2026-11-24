# circle_rank_map() の sigma のラプラス近似を MCMC と比べる ------------------------------------------------
# (1) sigma: 同時の最大値(MAP)・ラプラス近似の山・MCMC の平均・中央値・最頻値(log sigma の密度の山)
# (2) 置き直しで布置と半径がどれくらい動くか(MAP の解とラプラスの解の差，それぞれと MCMC の事後平均の差)
# データは 20261124_circle_rank_map_check.R と同じ 5 つ・2 次元と 3 次元。MCMC は adapt_delta 0.99

pacman::p_load(tidyverse, cmdstanr, posterior)
source("develop/circle_rank_map.R")
mod <- cmdstanr::cmdstan_model("develop/circle_rank.stan")

local({
  ex <- parse("develop/20261122_ellipse_rank_recovery_oi1987.R")
  for (e in ex) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("simulate_oi"))) {
      eval(e, envir = globalenv())
    }
  }
})

# MCMC の事後平均(点間距離・半径)と sigma の要約 ------------------------------------------------------
mcmc_summary <- function(dat) {
  od <- tempfile("fit_")
  dir.create(od)
  f <- mod$sample(data = dat[c("N", "P", "L", "row", "col", "K", "y")],
                  init = function(chain_id) .cr_init(dat$Z0, dat), chains = 4, parallel_chains = 4,
                  iter_warmup = 1000, iter_sampling = 2500, seed = 1, adapt_delta = 0.99,
                  refresh = 0, show_messages = FALSE, output_dir = od)
  N <- dat$N
  P <- dat$P
  dr <- posterior::as_draws_matrix(f$draws(c("Lx", "r", "sig")))
  # 各ドローの点間距離(平方和 N に正規化した布置)
  D <- t(apply(dr[, grepl("^Lx\\[", colnames(dr)), drop = FALSE], 1, function(v) {
    X <- rbind(0, matrix(v, N - 1, P))
    C <- sweep(X, 2, colMeans(X))
    c(dist(sqrt(N / sum(C^2)) * C))
  }))
  sig <- as.numeric(dr[, "sig"])
  dn <- density(log(sig))
  list(D = colMeans(D), r = colMeans(dr[, grepl("^r\\[", colnames(dr)), drop = FALSE]),
       sig = c(mean = mean(sig), median = median(sig), mode = exp(dn$x[which.max(dn$y)])))
}

check_one <- function(sim, P) {
  fm <- circle_rank_map(sim$d, ndim = P, sigma = "map")
  t0 <- Sys.time()
  fl <- circle_rank_map(sim$d, ndim = P, sigma = "laplace")
  sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  dat <- fl$dat
  mc <- mcmc_summary(dat)
  # 足場の並び(dat$ord)にそろえて比べる
  Dm <- c(dist(fm$conf[dat$ord, , drop = FALSE]))
  Dl <- c(dist(fl$conf[dat$ord, , drop = FALSE]))
  rm <- fm$radius[dat$ord]
  rl <- fl$radius[dat$ord]
  tibble(sig_map = fm$sigma, sig_lap = fl$sigma,
         sig_mean = mc$sig[["mean"]], sig_median = mc$sig[["median"]], sig_mode = mc$sig[["mode"]],
         dD_map_lap = max(abs(Dm - Dl)), dr_map_lap = max(abs(rm - rl)),
         rng_r_map = max(rm), rng_r_lap = max(rl),
         dD_map_mc = max(abs(Dm - mc$D)), dD_lap_mc = max(abs(Dl - mc$D)),
         dr_map_mc = max(abs(rm - mc$r)), dr_lap_mc = max(abs(rl - mc$r)),
         rng_r_mc = max(mc$r), sec_lap = sec)
}

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
    check_one(simulate_oi(N = N, R = R, E = E, seed = seed), P)
  })) %>%
  unnest(out)
