# フォンミーゼスモデル(AMISESCAL)のベイズ版を Stan で試す --------------------------------
# モデル: develop/vonmises_mds.stan(識別はコレスキー因子＝cholesky_factor_cov で縛る)
# 関数: develop/vonmises_functions.R(データ生成に .vm_parts を使う)

pacman::p_load(tidyverse, cmdstanr, posterior)
source("develop/vonmises_functions.R")

mod <- cmdstanr::cmdstan_model("develop/vonmises_mds.stan")
KAPPA_FIXED_MAX <- 6.5   # omega = 1/(2pi) で pi < 1 になる上限 6.5508 の手前

# 基準の3点を選んで先頭に並べる ---------------------------------------------------------
# 対称部の Torgerson 解で，重心から最も遠い点→そこから最も遠い点→2点を結ぶ直線から最も遠い点
choose_anchor <- function(G) {
  D <- G
  D[is.na(D)] <- t(D)[is.na(D)]
  D[is.na(D)] <- mean(D, na.rm = TRUE)
  diag(D) <- 0
  Z <- stats::cmdscale((D + t(D)) / 2, k = 2)
  a1 <- which.max(rowSums(sweep(Z, 2, colMeans(Z))^2))
  a2 <- which.max(rowSums(sweep(Z, 2, Z[a1, ])^2))
  u <- (Z[a2, ] - Z[a1, ]) / sqrt(sum((Z[a2, ] - Z[a1, ])^2))
  W <- sweep(Z, 2, Z[a1, ])
  a3 <- which.max(abs(W[, 1] * u[2] - W[, 2] * u[1]))
  c(a1, a2, a3, setdiff(seq_len(nrow(G)), c(a1, a2, a3)))
}

# Stan に渡すデータを作る ---------------------------------------------------------------
make_standata <- function(G, omega_fixed = 1, prior_x_sd = 1, prior_c = c(2, 4)) {
  ord <- choose_anchor(G)
  G <- G[ord, ord]
  diag(G) <- NA
  sc <- sqrt(mean(G^2, na.rm = TRUE))
  use <- !is.na(G)
  list(data = list(N = nrow(G), L = sum(use), row = row(G)[use], col = col(G)[use],
                   g = G[use] / sc, omega_fixed = omega_fixed,
                   kappa_max = if (omega_fixed == 1) KAPPA_FIXED_MAX else 20,
                   prior_x_sd = prior_x_sd, prior_c_a = prior_c[1], prior_c_b = prior_c[2]),
       order = ord, scale = sc, labels = rownames(G))
}

# 当てはめて要約する -------------------------------------------------------------------
fit_vm_stan <- function(G, omega_fixed = 1, seed = 1, nopt = 10, ...) {
  sd <- make_standata(G, omega_fixed = omega_fixed, ...)
  # 局所解が多いので，まず事後モードを乱数の初期値 nopt 本で探し，最良のものを全チェーンの初期値にする
  opts <- map(seq_len(nopt), function(k) {
    tryCatch(mod$optimize(data = sd$data, seed = seed * 100 + k, refresh = 0,
                          show_messages = FALSE, show_exceptions = FALSE, iter = 5000),
             error = function(e) NULL)
  }) %>% compact()
  lp <- map_dbl(opts, ~ .x$lp())
  best <- opts[[which.max(lp)]]
  init <- map(1:4, ~ list(Lx = unname(best$mle("Lx")) %>% matrix(sd$data$N - 1, 2),
                          m = map(seq_len(sd$data$N), function(r) {
                            v <- c(cos(unname(best$mle("mu"))[r]), sin(unname(best$mle("mu"))[r]))
                            v
                          }),
                          kappa = pmin(unname(best$mle("kappa")), sd$data$kappa_max * 0.99),
                          c_free = if (omega_fixed == 1) numeric(0) else pmin(pmax(unname(best$mle("c")), 0.01), 0.99),
                          sigma = unname(best$mle("sigma"))))
  fit <- mod$sample(data = sd$data, seed = seed, chains = 4, parallel_chains = 4,
                    iter_warmup = 1000, iter_sampling = 1000, refresh = 0,
                    show_messages = FALSE, init = init)
  list(fit = fit, sd = sd, opt_lp = sort(lp, decreasing = TRUE))
}

summarise_vm <- function(res) {
  s <- res$fit$summary(c("X", "mu", "kappa", "c", "sigma"),
                       "mean", "sd", ~ quantile(.x, c(0.025, 0.975)), "rhat", "ess_bulk")
  # mu は円周上の量なので平均は単位ベクトルの平均の向きで出す
  dr <- res$fit$draws("mu", format = "draws_matrix")
  mu_circ <- atan2(colMeans(sin(dr)), colMeans(cos(dr)))
  mu_R <- sqrt(colMeans(sin(dr))^2 + colMeans(cos(dr))^2)   # 平均合成ベクトル長(1 に近いほど確か)
  list(summary = s, mu_circ = mu_circ, mu_R = mu_R,
       diag = res$fit$diagnostic_summary(quiet = TRUE))
}

# 真値をコレスキーの枠(点1を原点，点2を +x 軸，点3を上半面)へ移す --------------------------------
to_chol_frame <- function(X, mu) {
  X <- sweep(X, 2, X[1, ])
  ang <- atan2(X[2, 2], X[2, 1])
  R <- matrix(c(cos(-ang), sin(-ang), -sin(-ang), cos(-ang)), 2)
  X <- X %*% t(R)
  mu <- mu - ang
  if (X[3, 2] < 0) {
    X[, 2] <- -X[, 2]
    mu <- -mu
  }
  list(X = X, mu = ((mu + pi) %% (2 * pi)) - pi)
}

# 1. パラメタリカバリ: omega 固定(2011 版)で作ったデータ --------------------------------------
sim_data <- function(n, sd_err, seed) {
  set.seed(seed)
  X <- matrix(stats::rnorm(2 * n), n, 2)
  mu <- stats::runif(n, -pi, pi)
  kappa <- stats::runif(n, 1, 6)
  peak <- exp(kappa) / (2 * pi * besselI(kappa, 0))
  cl <- list(n = n, i = row(diag(n))[row(diag(n)) != col(diag(n))],
             j = col(diag(n))[row(diag(n)) != col(diag(n))])
  # .vm_parts は kappa = KMAX * plogis(a) の母数で受ける
  par_true <- c(X, mu, stats::qlogis(kappa / KMAX), stats::qlogis(peak))
  xi <- .vm_parts(par_true, cl)$xi
  G <- matrix(NA_real_, n, n, dimnames = list(LETTERS[1:n], LETTERS[1:n]))
  G[cbind(cl$i, cl$j)] <- pmax(xi + stats::rnorm(length(xi), 0, sd_err), 0.01)
  list(G = G, X = X, mu = mu, kappa = kappa)
}

sim <- sim_data(n = 12, sd_err = 0.1, seed = 3)
res_sim <- fit_vm_stan(sim$G, omega_fixed = 1)
sm_sim <- summarise_vm(res_sim)
sm_sim$diag
round(res_sim$opt_lp, 1)
# 真値を同じ並べ替えとデータの単位，コレスキーの枠へ
ord <- res_sim$sd$order
tr <- to_chol_frame(sim$X[ord, ] / res_sim$sd$scale, sim$mu[ord])
X_post <- matrix(sm_sim$summary$mean[1:(2 * 12)], 12, 2)
recovery_sim <- tibble(
  label = res_sim$sd$labels,
  x1_true = tr$X[, 1], x1_post = X_post[, 1],
  x2_true = tr$X[, 2], x2_post = X_post[, 2],
  mu_true_deg = tr$mu * 180 / pi, mu_post_deg = sm_sim$mu_circ * 180 / pi, mu_R = sm_sim$mu_R,
  kappa_true = sim$kappa[ord],
  kappa_post = sm_sim$summary %>% filter(str_starts(variable, "kappa")) %>% pull(mean)
)
recovery_sim %>% mutate(across(where(is.numeric), ~ round(.x, 2)))
sm_sim$summary %>% summarise(max_rhat = max(rhat, na.rm = TRUE), min_ess = min(ess_bulk, na.rm = TRUE))

# 2. 千野 (1997, p.13) のソシオメトリックデータ ---------------------------------------------
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

res_fix <- fit_vm_stan(chino, omega_fixed = 1)
res_free <- fit_vm_stan(chino, omega_fixed = 0)
sm_fix <- summarise_vm(res_fix)
sm_free <- summarise_vm(res_free)
list(fixed = sm_fix$diag, free = sm_free$diag)
list(fixed = round(res_fix$opt_lp, 1), free = round(res_free$opt_lp, 1))
map_dfr(list(fixed = sm_fix, free = sm_free),
        ~ .x$summary %>% summarise(max_rhat = max(rhat, na.rm = TRUE), min_ess = min(ess_bulk, na.rm = TRUE)),
        .id = "omega")
# 各人の平均方向と確かさ(R)・集中度・山の高さ
chino_table <- function(res, sm) {
  s <- sm$summary
  tibble(label = res$sd$labels,
         mu_deg = round(sm$mu_circ * 180 / pi), mu_R = round(sm$mu_R, 2),
         kappa = round(s %>% filter(str_starts(variable, "kappa")) %>% pull(mean), 2),
         c = round(s %>% filter(str_starts(variable, stringr::fixed("c["))) %>% pull(mean), 2))
}
chino_table(res_fix, sm_fix)
chino_table(res_free, sm_free)
# LOO で omega 固定と自由を比べる
loo_fix <- res_fix$fit$loo()
loo_free <- res_free$fit$loo()
loo::loo_compare(list(fixed = loo_fix, free = loo_free))

# 3. 図: 事後平均の布置に，事後の各ドローの mu の向きを細い線で重ねる ---------------------------------
plot_vm_post <- function(res, sm, main = "") {
  s <- sm$summary
  n <- res$sd$data$N
  X <- matrix(s$mean[1:(2 * n)], n, 2)
  dr <- res$fit$draws("mu", format = "draws_matrix")
  rho <- 0.12 * diff(range(X))
  lim <- range(X) + c(-1.6, 1.6) * rho
  plot(X, asp = 1, xlim = lim, ylim = lim, type = "n", xlab = "Dim 1", ylab = "Dim 2", main = main)
  abline(h = 0, v = 0, col = "grey80", lty = 3)
  take <- seq(1, nrow(dr), length.out = 200)
  for (r in seq_len(n)) {
    segments(X[r, 1], X[r, 2], X[r, 1] + rho * cos(dr[take, r]), X[r, 2] + rho * sin(dr[take, r]),
             col = grDevices::adjustcolor("firebrick", 0.05))
  }
  points(X, pch = 21, bg = "white", cex = 2.2)
  text(X, labels = res$sd$labels, cex = 0.8)
}
out_png <- file.path(Sys.getenv("VM_OUT", tempdir()), "vonmises_stan.png")
png(out_png, width = 1500, height = 520, res = 110)
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
plot_vm_post(res_sim, sm_sim, "simulation (omega fixed)")
plot_vm_post(res_fix, sm_fix, "Chino: omega fixed")
plot_vm_post(res_free, sm_free, "Chino: omega free")
dev.off()
out_png
