# circle_rank_map()(base R の事後確率最大値・develop/circle_rank_map.R)を Stan 版と照合する ----------------------
# (1) 対数事後密度が Stan の optimize の lp__ と同じ値か(Stan の最大値の点で R の関数を評価する)
# (2) 解析的な勾配が数値微分と合うか
# (3) R の最大化の答えが Stan の optimize と同じ山か(対数事後密度・布置・半径・sigma)と，時間
# データは 20261123_circle_rank_scaffold.R と同じ 5 つ・2 次元と 3 次元

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

# Stan の最大値の点を R の母数の形にする ----------------------------------------------------------
stan_par <- function(o, dat) {
  m <- o$mle()
  list(Lx = matrix(m[grepl("^Lx\\[", names(m))], dat$N - 1, dat$P),
       a_raw = unname(m[grepl("^a_raw\\[", names(m))]), s_a = m[["s_a"]], sig = m[["sig"]],
       g = unname(m[grepl("^g\\[", names(m))]))
}

# 数値微分(中心差分) ---------------------------------------------------------------------------
num_grad <- function(theta, dat, h = 1e-6) {
  vapply(seq_along(theta), function(i) {
    e <- replace(numeric(length(theta)), i, h)
    (.cr_logpost(.cr_unpack(theta + e, dat), dat) - .cr_logpost(.cr_unpack(theta - e, dat), dat)) / (2 * h)
  }, numeric(1))
}

check_one <- function(sim, P) {
  # R の最大化
  t0 <- Sys.time()
  fr <- circle_rank_map(sim$d, ndim = P, nrandom = 10)
  sec_r <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  dat <- fr$dat
  # Stan の最大化(同じデータ・同じ並び・Torgerson 解から)
  sdat <- dat[c("N", "P", "L", "row", "col", "K", "y")]
  od <- tempfile("fit_")
  dir.create(od)
  t0 <- Sys.time()
  o <- mod$optimize(data = sdat, init = list(.cr_init(dat$Z0, dat)), seed = 1, jacobian = FALSE,
                    iter = 20000, refresh = 0, show_messages = FALSE, output_dir = od)
  sec_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ps <- stan_par(o, dat)
  # (1) 同じ点での対数事後密度
  lp_r_at_stan <- .cr_logpost(ps, dat)
  # (2) 勾配(Stan の最大値の点と，R の初期値の点)
  th_s <- .cr_pack(ps, dat)
  th_0 <- .cr_pack(.cr_init(dat$Z0, dat), dat)
  gerr <- max(abs(.cr_grad_theta(th_0, dat) - num_grad(th_0, dat)))
  gmax_stan <- max(abs(.cr_grad_theta(th_s, dat)))   # Stan の最大値の点で勾配が 0 に近いか
  # (3) 答えの比べ合い(点間距離・半径・sigma)
  N <- dat$N
  X <- rbind(0, ps$Lx)
  C <- sweep(X, 2, colMeans(X))
  tau_s <- sqrt(N / sum(C^2)) * C
  D_s <- c(dist(tau_s))
  D_r <- c(dist(fr$conf[dat$ord, , drop = FALSE]))
  a_s <- ps$s_a * ps$a_raw
  r_s <- a_s - min(a_s)
  tibble(lp_stan = o$lp(), lp_r_at_stan = lp_r_at_stan, lp_r = fr$lp,
         grad_err = gerr, grad_at_stan = gmax_stan,
         max_dD = max(abs(D_r - D_s)), max_dr = max(abs(fr$radius[dat$ord] - r_s)),
         sig_r = fr$sigma, sig_stan = ps$sig,
         n_same = fr$n_same, n_start = length(fr$lp_starts), conv = fr$convergence,
         sec_r = sec_r, sec_stan = sec_s)
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
