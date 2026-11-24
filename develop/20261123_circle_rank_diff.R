# 円モデルの半径を「差だけ」の母数にした版(circle_rank.stan・2026-11-23)を前の版と比べる -------------------
# mu = d - r_j + r_k は全員の半径に同じ数を足しても変わらない。前の版は r = exp(c)・c = m_c + s_c * c_raw の
# 対数正規で半径の水準を事前分布だけで決めていて，m_c と s_c の事後分布が尾根になった(相関 -0.92)。
# 新しい版は和が 0 の a_j ~ N(0, s_a^2) を推定し，OI と同じく r_j = a_j - min(a) と最小を 0 にずらして報告する。
# 同じ 5 つのデータ・1〜3 次元で，(1) 尾根が消えて有効標本数が増えるか，(2) 発散が減るか，
# (3) 布置・半径・sigma の回復と周辺尤度の次元の選び方が変わるかを見る
# 前の版の Stan は git の 5860065 の develop/circle_rank.stan(ここでは一時フォルダに書き出して使う)

pacman::p_load(tidyverse, cmdstanr, posterior, bridgesampling)

# 前の版と新しい版の Stan を用意する ---------------------------------------------------------------
old_stan <- file.path(tempdir(), "circle_rank_lognormal.stan")
writeLines(system2("git", c("show", "5860065:develop/circle_rank.stan"), stdout = TRUE), old_stan)
mods <- list(
  lognormal = cmdstanr::cmdstan_model(old_stan, compile_model_methods = TRUE, force_recompile = TRUE),
  diff = cmdstanr::cmdstan_model("develop/circle_rank.stan", compile_model_methods = TRUE,
                                 force_recompile = TRUE)
)

# simulate_oi() と，足場のスクリプトの pick_scaffold()・to_scaffold()・prep_circle() を読み込む -------------
load_defs <- function(file, names) {
  ex <- parse(file)
  for (e in ex) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]]) %in% names) {
      eval(e, envir = globalenv())
    }
  }
}
load_defs("develop/20261122_ellipse_rank_recovery_oi1987.R", "simulate_oi")
load_defs("develop/20261123_circle_rank_scaffold.R", c("pick_scaffold", "to_scaffold", "prep_circle"))

# 当てはめて比べる量を出す ---------------------------------------------------------------------------
fit_one <- function(sim, P, model) {
  pr <- prep_circle(sim$d, P)
  N <- pr$dat$N
  init <- function(chain_id) {
    i0 <- pr$make_init(chain_id)
    if (model == "diff") {
      i0[c("c_raw", "m_c", "s_c")] <- NULL
      i0$a_raw <- rep(0, N)
      i0$s_a <- 0.3
    }
    i0
  }
  D_true <- sim$D[pr$ord, pr$ord]
  od <- tempfile("fit_")
  dir.create(od)
  t0 <- Sys.time()
  f <- mods[[model]]$sample(data = pr$dat, init = init, chains = 4, parallel_chains = 4,
                            iter_warmup = 1000, iter_sampling = 2500, seed = 1, refresh = 0,
                            show_messages = FALSE, output_dir = od, sig_figs = 18)
  sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  # 有効標本数と Rhat は，推定する母数(座標・半径・sig・閾値)と lp__ について
  pars <- if (model == "diff") c("lp__", "tau", "a", "s_a", "sig", "g") else
    c("lp__", "tau", "c", "m_c", "s_c", "sig", "g")
  sm <- f$summary(pars, "rhat", "ess_bulk", "ess_tail")
  rad_par <- if (model == "diff") c("a", "s_a") else c("c", "m_c", "s_c")
  sm_r <- f$summary(rad_par, "ess_bulk")
  tau <- unclass(as.matrix(f$draws("tau", format = "matrix")))
  Dhat <- colMeans(t(apply(tau, 1, function(v) c(dist(matrix(v, N, P))))))
  # 半径: 新しい版は r = a - min(a)，前の版は exp(c)。水準は事後分布の事前分布任せなので，比べるのは差の回復
  r_hat <- if (model == "diff") colMeans(unclass(as.matrix(f$draws("r", format = "matrix")))) else
    colMeans(exp(unclass(as.matrix(f$draws("c", format = "matrix")))))
  # 尾根の目安: 前の版は m_c と s_c，新しい版は s_a と sig の事後相関
  ridge <- if (model == "diff") cor(c(f$draws("s_a")), c(f$draws("sig"))) else
    cor(c(f$draws("m_c")), c(f$draws("s_c")))
  f$init_model_methods()
  U <- unclass(as.matrix(f$unconstrain_draws(format = "draws_matrix")))
  colnames(U) <- paste0("u", seq_len(ncol(U)))
  lb <- setNames(rep(-Inf, ncol(U)), colnames(U))
  ub <- setNames(rep(Inf, ncol(U)), colnames(U))
  lp <- function(x, data) {
    tryCatch(f$log_prob(unconstrained_variables = x, jacobian = TRUE), error = function(e) -Inf)
  }
  b <- bridge_sampler(samples = U, log_posterior = lp, data = NULL, lb = lb, ub = ub,
                      method = "warp3", repetitions = 3, silent = TRUE)
  tibble(logml = mean(b$logml),
         logml_sd = sd(b$logml),
         divergent = sum(f$diagnostic_summary(quiet = TRUE)$num_divergent),
         rhat_max = max(sm$rhat),
         ess_min = min(sm$ess_bulk),
         ess_tail_min = min(sm$ess_tail),
         ess_radius_min = min(sm_r$ess_bulk),
         ridge_cor = ridge,
         sigma = mean(f$draws("sig")),
         cor_dist = cor(Dhat, c(D_true[lower.tri(D_true)])),
         cor_r = cor(r_hat, sim$r[pr$ord]),
         sec = sec)
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
res <- dsets %>%
  mutate(out = pmap(list(N, R, E, seed), function(N, R, E, seed) {
    sim <- simulate_oi(N = N, R = R, E = E, seed = seed)
    expand_grid(model = c("lognormal", "diff"), P = 1:3) %>%
      mutate(fit = map2(model, P, function(m, P) fit_one(sim, P, m))) %>%
      unnest(fit)
  })) %>%
  unnest(out)
saveRDS(res, "develop/20261123_circle_rank_diff.rds")
