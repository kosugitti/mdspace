# ellipse_rank_map()(base R の事後確率最大値・develop/ellipse_rank_map.R)を Stan 版(ellipse_rank_std.stan)と照合する ---------
# (1) 対数事後密度: Stan の optimize の最大値の点で R の関数を評価し，lp__ と同じ値か
# (2) 勾配: Stan の最大値の点で R の勾配が 0 に近いか
# (3) 最大化の答え: R(Torgerson 解＋乱数 10 本)と Stan の optimize(Torgerson 解から)の対数事後密度・布置・半軸・sigma
# データ: 職業データ(書籍表1.19・類似度)と，模擬データ 2 つ(20261122_ellipse_rank_recovery_oi1987.R の simulate_oi)
# 下位モデル 4 つ・2 次元

pacman::p_load(tidyverse, cmdstanr, posterior)
source("develop/ellipse_rank_map.R")
mod <- cmdstanr::cmdstan_model("develop/ellipse_rank_std.stan")

local({
  ex <- parse("develop/20261122_ellipse_rank_recovery_oi1987.R")
  for (e in ex) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("simulate_oi"))) {
      eval(e, envir = globalenv())
    }
  }
})

occ <- as.matrix(read.csv("tests/testthat/fixtures/okada_imaizumi_2024_table1_19.csv",
                          row.names = 1, check.names = FALSE))

# Stan の最大値の点を R の母数の形にする ----------------------------------------------------------
stan_par <- function(o, dat) {
  m <- o$mle()
  N <- dat$N
  P <- dat$P
  get <- function(nm) unname(m[grepl(paste0("^", nm, "\\["), names(m))])
  list(yr = matrix(get("yr"), N - 1, P), v = get("v"),
       c_raw = if (dat$use_size) get("c_raw") else rep(0, N),
       et = if (dat$use_shape) matrix(get("et"), N, P, byrow = FALSE) else matrix(0, N, P),
       m_c = m[["m_c"]], s_c = if (dat$use_size) m[["s_c"]] else 0,
       s_e = if (dat$use_shape) m[["s_e"]] else 0, sig = m[["sig"]], g = get("g"))
}

check_one <- function(x, type, label, model, P = 2) {
  t0 <- Sys.time()
  fr <- ellipse_rank_map(x, ndim = P, type = type, model = model, nrandom = 10)
  sec_r <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  dat <- fr$dat
  sdat <- c(dat[c("N", "P", "L", "row", "col", "K", "y", "anchor")],
            list(use_size = as.integer(dat$use_size), use_shape = as.integer(dat$use_shape)))
  ini <- .er_init(dat$Z0, dat)
  # Stan では使わない母数も事前分布で動くので，初期値は小さい正の値に
  sini <- list(yr = ini$yr, v = ini$v, c_raw = ini$c_raw, et = lapply(seq_len(dat$N), function(j) ini$et[j, ]),
               m_c = ini$m_c, s_c = max(ini$s_c, 0.1), s_e = max(ini$s_e, 0.1), sig = ini$sig, g = ini$g)
  od <- tempfile("fit_")
  dir.create(od)
  t0 <- Sys.time()
  o <- mod$optimize(data = sdat, init = list(sini), seed = 1, jacobian = FALSE, iter = 50000,
                    refresh = 0, show_messages = FALSE, output_dir = od)
  sec_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ps <- stan_par(o, dat)
  lp_stan <- o$lp()
  # Stan は使わない母数(s_c・s_e)を 0 の手前で止めるので，R の側では使う母数だけで評価する
  ps_r <- ps
  if (!dat$use_size) { ps_r$c_raw <- rep(0, dat$N); ps_r$s_c <- 0 }
  if (!dat$use_shape) { ps_r$et <- matrix(0, dat$N, dat$P); ps_r$s_e <- 0 }
  lp_r_at_stan <- .er_logpost(ps, dat)
  gmax_stan <- max(abs(.er_grad_theta(.er_pack(ps_r, dat), dat)))
  # 答え
  rs <- .er_result(.er_pack(ps_r, dat), dat, rownames(fr$conf))
  tibble(data = label, model = model,
         lp_stan = lp_stan, lp_diff = lp_r_at_stan - lp_stan, lp_r = fr$lp,
         grad_at_stan = gmax_stan,
         cor_dist = cor(c(dist(fr$conf)), c(dist(rs$conf))),
         conf_maxdiff = max(abs(fr$conf - rs$conf)),
         axis_maxdiff = max(abs(log(fr$semiaxis) - log(rs$semiaxis))),
         sigma_r = fr$sigma, sigma_stan = ps$sig,
         neg_mu = fr$negative_mu, sec_r = sec_r, sec_stan = sec_s, best_start = which.max(fr$lp_starts))
}

sims <- list(
  list(x = occ, type = "similarity", label = "職業(表1.19)"),
  list(x = simulate_oi(N = 12, R = 0.3, E = 0.5, seed = 2001)$d, type = "dissimilarity", label = "N12・E0.5"),
  list(x = simulate_oi(N = 10, R = 0.3, E = 0.3, seed = 1007)$d, type = "dissimilarity", label = "N10・E0.3")
)
models <- c("ellipse", "circle", "equal_area", "symmetric")
res <- map_dfr(sims, function(s) map_dfr(models, function(m) check_one(s$x, s$type, s$label, m)))
saveRDS(res, file.path(tempdir(), "er_check.rds")); options(width = 250); res %>% mutate(across(where(is.double), ~ signif(.x, 4))) %>% as.data.frame() %>% print()
