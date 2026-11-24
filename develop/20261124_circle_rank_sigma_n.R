# 対象の数 N を増やすと，sigma の MAP と MCMC のずれは縮むか ------------------------------------------------
# 予想: セルは N(N-1) で増え，布置・半径の母数は N で増えるので，母数 1 つあたりの情報が増えてずれは縮む
# N = 12, 20, 30, 50・E = 0.5・2 次元・各 2 データ。10 段階に粗くしたもの(閾値の数は一定)と粗くしないもの
# MCMC は adapt_delta 0.95・warmup 1000・sampling 1500(N = 50 の時間を抑えるため)

pacman::p_load(tidyverse, cmdstanr, posterior)
source("develop/circle_rank_map.R")
mod <- cmdstanr::cmdstan_model("develop/circle_rank.stan")

local({
  for (f in c("develop/20261122_ellipse_rank_recovery_oi1987.R", "develop/20261124_circle_rank_sigma_map.R")) {
    for (e in parse(f)) {
      if (is.call(e) && identical(e[[1]], as.name("<-")) &&
          as.character(e[[2]])[1] %in% c("simulate_oi", "coarsen")) {
        eval(e, envir = globalenv())
      }
    }
  }
})

run_one <- function(N, seed, K) {
  sim <- simulate_oi(N = N, R = 0.3, E = 0.5, seed = seed)
  d <- coarsen(sim$d, K)
  t0 <- Sys.time()
  fm <- circle_rank_map(d, ndim = 2, nrandom = 5, sigma = "map")
  sec_map <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  dat <- fm$dat
  od <- tempfile("fit_")
  dir.create(od)
  t0 <- Sys.time()
  f <- mod$sample(data = dat[c("N", "P", "L", "row", "col", "K", "y")],
                  init = function(chain_id) .cr_init(dat$Z0, dat), chains = 4, parallel_chains = 4,
                  iter_warmup = 1000, iter_sampling = 1500, seed = 1, adapt_delta = 0.95,
                  refresh = 0, show_messages = FALSE, output_dir = od)
  sec_mc <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  sig <- c(f$draws("sig"))
  r_mc <- colMeans(posterior::as_draws_matrix(f$draws("r")))
  tibble(K_used = dat$K, n_par = 2 * (N - 1) - 1 + (N - 1), L = dat$L,
         sig_map = fm$sigma, sig_mean = mean(sig), sig_median = median(sig),
         ratio = fm$sigma / mean(sig),
         r_cor = cor(fm$radius[dat$ord], r_mc), r_max_ratio = max(fm$radius) / max(r_mc),
         div = sum(f$sampler_diagnostics(format = "draws_df")$divergent__),
         sec_map = sec_map, sec_mc = sec_mc)
}

res <- expand_grid(N = c(12, 20, 30, 50), seed = c(3001, 3002), K = c(10, Inf)) %>%
  mutate(out = pmap(list(N, seed, K), run_one)) %>%
  unnest(out)
