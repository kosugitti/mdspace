# vonmises_mds()(base R の MAP)を Stan の事後モードと照合する ----------------------------
# Stan: develop/vonmises_mds.stan を $optimize()(ヤコビアンなし)で初期値20本，最良の lp__ を取る
# 比べるもの: ストレス・sigma・mu・kappa・c(布置は枠が違うので距離行列で比べる)

pacman::p_load(tidyverse, cmdstanr)
devtools::load_all(quiet = TRUE)

mod <- cmdstanr::cmdstan_model("develop/vonmises_mds.stan")

chino <- matrix(c(
  NA, 5, 4, 5, 3, 3, 2, 4, 1, 1,
  6, NA, 4, 2, 1, 2, 3, 3, 4, 3,
  4, 4, NA, 3, 3, 4, 4, 5, 4, 3,
  4, 1, 2, NA, 1, 1, 4, 2, 4, 3,
  7, 1, 2, 1, NA, 1, 2, 2, 2, 3,
  4, 3, 4, 2, 3, NA, 4, 4, 4, 4,
  4, 3, 4, 4, 5, 5, NA, 2, 4, 2,
  6, 4, 4, 4, 3, 4, 3, NA, 4, 4,
  2, 3, 3, 3, 3, 2, 3, 3, NA, 2,
  4, 4, 4, 5, 4, 4, 4, 4, 4, NA), 10, byrow = TRUE,
  dimnames = list(LETTERS[1:10], LETTERS[1:10]))

# Stan の事後モード(並べ替えは vonmises_mds と同じ基準点で) ---------------------------------
stan_map <- function(G, omega_fixed, nopt = 20) {
  ord <- mdspace:::.vm_anchor({g <- G; diag(g) <- NA; g})
  Gs <- G[ord, ord]
  diag(Gs) <- NA
  sc <- sqrt(mean(Gs^2, na.rm = TRUE))
  use <- !is.na(Gs)
  dat <- list(N = nrow(G), L = sum(use), row = row(Gs)[use], col = col(Gs)[use],
              g = Gs[use] / sc, omega_fixed = omega_fixed,
              kappa_max = if (omega_fixed == 1) 6.5 else 20,
              prior_x_sd = 1, prior_c_a = 2, prior_c_b = 4)
  opts <- map(seq_len(nopt), function(k) {
    tryCatch(mod$optimize(data = dat, seed = k, refresh = 0, iter = 10000,
                          show_messages = FALSE, show_exceptions = FALSE),
             error = function(e) NULL)
  }) %>% compact()
  best <- opts[[which.max(map_dbl(opts, ~ .x$lp()))]]
  X <- matrix(unname(best$mle("X")), nrow(G), 2)[order(ord), ]
  xi <- unname(best$mle("xi"))
  list(D = as.matrix(dist(X)), mu_rel = NULL,
       stress = sqrt(sum((dat$g - xi)^2) / sum(dat$g^2)),
       sigma = unname(best$mle("sigma")),
       kappa = unname(best$mle("kappa"))[order(ord)],
       c = unname(best$mle("c"))[order(ord)],
       mu = unname(best$mle("mu"))[order(ord)], X = X)
}

# 角度の差(回転・鏡映をそろえてから): mu を布置と同じ枠に移して比べる ---------------------------
angle_diff <- function(X1, mu1, X2, mu2) {
  # X2 を X1 に Procrustes(回転・鏡映)で合わせ，mu2 も同じ変換で動かす
  A <- sweep(X1, 2, colMeans(X1))
  B <- sweep(X2, 2, colMeans(X2))
  s <- svd(t(B) %*% A)
  R <- s$u %*% t(s$v)
  v <- cbind(cos(mu2), sin(mu2)) %*% R
  mu2r <- atan2(v[, 2], v[, 1])
  abs(((mu2r - mu1 + pi) %% (2 * pi)) - pi) * 180 / pi
}

compare_one <- function(G, omega) {
  r <- vonmises_mds(G, omega = omega, nrandom = 20)
  s <- stan_map(G, omega_fixed = as.integer(omega == "fixed"))
  tibble(omega = omega,
         stress_R = r$stress, stress_Stan = s$stress,
         sigma_R = r$sigma, sigma_Stan = s$sigma,
         maxdiff_dist = max(abs(as.matrix(dist(r$conf)) - s$D)),
         maxdiff_kappa = max(abs(r$kappa - s$kappa)),
         maxdiff_c = max(abs(r$c - s$c)),
         maxdiff_mu_deg = max(angle_diff(r$conf, r$mu, s$X, s$mu)))
}

check <- map_dfr(c("fixed", "free"), ~ compare_one(chino, .x))
check %>% mutate(across(where(is.numeric), ~ signif(.x, 4)))

# 同じ目的関数(R の負の対数事後密度)で両方の解を評価する --------------------------------------
eval_both <- function(G, omega) {
  r <- vonmises_mds(G, omega = omega, nrandom = 20)
  s <- stan_map(G, omega_fixed = as.integer(omega == "fixed"))
  ord <- mdspace:::.vm_anchor({g <- G; diag(g) <- NA; g})
  Gs <- G[ord, ord]
  diag(Gs) <- NA
  Gs <- Gs / sqrt(mean(Gs^2, na.rm = TRUE))
  use <- !is.na(Gs)
  free <- omega == "free"
  cells <- list(n = nrow(G), i = row(Gs)[use], j = col(Gs)[use], g = Gs[use], free = free,
                kmax = if (free) 20 else 6.5, prior = r$prior)
  ps <- mdspace:::.vm_pack(s$X[ord, ], s$mu[ord], s$kappa[ord], s$c[ord], s$sigma, cells)
  pr <- mdspace:::.vm_pack(r$conf[ord, ], r$mu[ord], r$kappa[ord], r$c[ord], r$sigma, cells)
  # Stan の解から R の optim を続けたらどこへ行くか
  cont <- optim(ps, mdspace:::.vm_nlp, mdspace:::.vm_gr, cells = cells, method = "BFGS",
                control = list(maxit = 5000, reltol = 1e-12))
  tibble(omega = omega, nlp_R = mdspace:::.vm_nlp(pr, cells), nlp_Stan = mdspace:::.vm_nlp(ps, cells),
         nlp_Stan_then_R = cont$value, R_starts_best5 = paste(round(sort(-r$starts)[1:5], 2), collapse = " "))
}
map_dfr(c("fixed", "free"), ~ eval_both(chino, .x))
