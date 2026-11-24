# 布置・半径を「近似の事後平均」で求めると MCMC に近づくか(INLA と同じ考え方) ------------------------------
# sigma の格子の各点で，他の母数(制約のない空間)を正規分布 N(theta_hat(sigma), H^-1) で近似し，
# そこから乱数を引いて点間距離と半径を計算する。sigma の周辺事後密度(ラプラス近似)の重みで平均する。
# あわせて，見立ての切り分けのため，各 sigma での正規近似の「中心」の半径も sigma の重みで平均する
# 格子の点でヘッセ行列が正定値にならないときは，その点を除いて平均する(除いた数は n_drop)
# (中心で MCMC に近ければ非対称の効果は小さい，中心でも膨らむなら最高点で代表させること自体の問題)

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

# 制約のない母数から，平方和 N の布置の点間距離と半径(最小を 0) -------------------------------------------
dist_radius <- function(theta, dat) {
  par <- .cr_unpack(theta, dat)
  X <- rbind(0, par$Lx)
  C <- sweep(X, 2, colMeans(X))
  a <- par$s_a * par$a_raw
  list(D = c(dist(sqrt(dat$N / sum(C^2)) * C)), r = a - min(a))
}

# MCMC の事後平均 ----------------------------------------------------------------------------------
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
  D <- t(apply(dr[, grepl("^Lx\\[", colnames(dr)), drop = FALSE], 1, function(v) {
    X <- rbind(0, matrix(v, N - 1, P))
    C <- sweep(X, 2, colMeans(X))
    c(dist(sqrt(N / sum(C^2)) * C))
  }))
  list(D = colMeans(D), r = colMeans(dr[, grepl("^r\\[", colnames(dr)), drop = FALSE]),
       sig = mean(dr[, "sig"]))
}

# 近似の事後平均 -----------------------------------------------------------------------------------
inla_mean <- function(fit, ngrid = 15, half = 0.5, ndraw = 1000, seed = 1) {
  dat <- fit$dat
  isig <- .cr_sig_index(dat)
  # 同時の最大値から始める(fit は sigma = "map" で当てはめたもの)
  cur <- fit$theta[-isig]
  ls_hat <- log(circle_rank_map_lap_sigma(fit))
  grid <- ls_hat + seq(-half, half, length.out = ngrid)
  set.seed(seed)
  pts <- lapply(grid, function(ls) {
    f <- .cr_fix_sig(ls, dat, jac = TRUE)
    r <- stats::optim(cur, f$fn, f$gr, method = "BFGS", control = list(maxit = 5000, reltol = 1e-10))
    cur <<- r$par
    H <- stats::optimHess(r$par, f$fn, f$gr)
    H <- (H + t(H)) / 2
    ld <- determinant(H, logarithm = TRUE)
    lm <- -r$value - 0.5 * as.numeric(ld$modulus)
    # 中心の距離・半径と，正規近似から引いた乱数の平均
    ctr <- dist_radius(f$full(r$par), dat)
    R <- tryCatch(chol(H), error = function(e) NULL)
    if (is.null(R)) return(NULL)                       # 正定値でない点は除く(数を記録)
    Z <- matrix(stats::rnorm(ndraw * length(cur)), ndraw)
    V <- sweep(t(backsolve(R, t(Z))), 2, r$par, `+`)
    dr <- lapply(seq_len(ndraw), function(i) dist_radius(f$full(V[i, ]), dat))
    list(lm = lm, D_c = ctr$D, r_c = ctr$r,
         D_m = rowMeans(sapply(dr, `[[`, "D")), r_m = rowMeans(sapply(dr, `[[`, "r")))
  })
  ok <- !vapply(pts, is.null, logical(1))
  n_drop <- sum(!ok)
  pts <- pts[ok]
  grid <- grid[ok]
  ngrid <- length(grid)
  lm <- sapply(pts, `[[`, "lm")
  w <- exp(lm - max(lm))
  w <- w / sum(w)
  wavg <- function(nm) Reduce(`+`, Map(function(p, wi) wi * p[[nm]], pts, w))
  list(D_c = wavg("D_c"), r_c = wavg("r_c"), D_m = wavg("D_m"), r_m = wavg("r_m"),
       sig = sum(w * exp(grid)), w_end = max(w[c(1, ngrid)]), n_drop = n_drop)
}
circle_rank_map_lap_sigma <- function(fit) {
  exp(.cr_sigma_laplace(fit$theta, fit$dat)$log_sigma)
}

check_one <- function(sim, P) {
  fm <- circle_rank_map(sim$d, ndim = P, sigma = "map")
  dat <- fm$dat
  dm <- dist_radius(fm$theta, dat)
  t0 <- Sys.time()
  il <- inla_mean(fm)
  sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  mc <- mcmc_summary(dat)
  tibble(sig_map = fm$sigma, sig_inla = il$sig, sig_mc = mc$sig, w_end = il$w_end, n_drop = il$n_drop,
         rmax_map = max(dm$r), rmax_ctr = max(il$r_c), rmax_inla = max(il$r_m), rmax_mc = max(mc$r),
         dr_map = max(abs(dm$r - mc$r)), dr_ctr = max(abs(il$r_c - mc$r)), dr_inla = max(abs(il$r_m - mc$r)),
         dD_map = max(abs(dm$D - mc$D)), dD_ctr = max(abs(il$D_c - mc$D)), dD_inla = max(abs(il$D_m - mc$D)),
         sec_inla = sec)
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
