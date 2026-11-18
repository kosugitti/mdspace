# フォンミーゼスモデル(AMISESCAL; 荘島 2011, 2012)の試作: optim 版の確認 --------------
# 関数は develop/vonmises_functions.R

pacman::p_load(tidyverse)
source("develop/vonmises_functions.R")

# 1. 勾配の確認(数値微分と比べる) --------------------------------------------------------
set.seed(11)
n <- 6
Gt <- matrix(stats::runif(n * n, 1, 7), n, n)
diag(Gt) <- NA
Gs <- Gt / sqrt(mean(Gt^2, na.rm = TRUE))
use <- !is.na(Gs)
cells <- list(n = n, i = row(Gs)[use], j = col(Gs)[use], g = Gs[use],
              love = pmax(t(Gs)[use] - Gs[use], 0))
par0 <- .vm_start(Gs, cells, random = TRUE)
num_grad <- function(f, x, ...) {
  vapply(seq_along(x), function(k) {
    h <- 1e-6
    e <- replace(numeric(length(x)), k, h)
    (f(x + e, ...) - f(x - e, ...)) / (2 * h)
  }, numeric(1))
}
grad_check <- tibble(
  lambda = c(0, 0.5),
  max_abs_diff = c(max(abs(.vm_grad(par0, cells, 0) - num_grad(.vm_obj, par0, cells = cells, lambda = 0))),
                   max(abs(.vm_grad(par0, cells, 0.5) - num_grad(.vm_obj, par0, cells = cells, lambda = 0.5))))
)
grad_check

# 2. パラメタリカバリ(lambda = 0 の 2011 版) ----------------------------------------------
# 真値から xi を作り，誤差を足して当てはめる。回転と鏡映は Procrustes で合わせる
recover_once <- function(n, sd_err, seed) {
  set.seed(seed)
  X <- matrix(stats::rnorm(2 * n), n, 2)
  mu <- stats::runif(n, -pi, pi)
  kappa <- stats::runif(n, 1, 6)
  peak <- stats::runif(n, 0.1, 0.7)
  par_true <- c(X, mu, stats::qlogis(kappa / KMAX), stats::qlogis(peak))
  cl <- list(n = n, i = row(diag(n))[row(diag(n)) != col(diag(n))],
             j = col(diag(n))[row(diag(n)) != col(diag(n))])
  xi <- .vm_parts(par_true, cl)$xi
  G <- matrix(NA_real_, n, n)
  G[cbind(cl$i, cl$j)] <- pmax(xi + stats::rnorm(length(xi), 0, sd_err), 0.01)
  fit <- vonmises_mds_try(G, nrandom = 10, seed = seed)
  # 真の布置をデータの単位(fit$scale で割る)にして Procrustes(回転・鏡映)で合わせる
  Xt <- sweep(X, 2, colMeans(X)) / fit$scale
  s <- svd(t(fit$conf) %*% Xt)
  R <- s$u %*% t(s$v)
  Xf <- fit$conf %*% R
  # mu も同じ回転・鏡映で動かす
  muv <- cbind(cos(fit$mu), sin(fit$mu)) %*% R
  mu_f <- atan2(muv[, 2], muv[, 1])
  dmu <- abs(((mu_f - mu + pi) %% (2 * pi)) - pi)
  tibble(n = n, sd_err = sd_err, seed = seed, stress = fit$stress,
         rmse_X = sqrt(mean((Xf - Xt)^2)),
         med_dmu_deg = stats::median(dmu) * 180 / pi,
         cor_kappa = stats::cor(log(fit$kappa), log(kappa)),
         cor_peak = stats::cor(fit$peak, peak))
}
recovery <- expand_grid(n = c(10, 20), sd_err = c(0, 0.1), seed = 1:3) %>%
  pmap_dfr(recover_once)
recovery %>%
  group_by(n, sd_err) %>%
  summarise(across(stress:cor_peak, ~ round(mean(.x), 3)), .groups = "drop")

# 3. スライドの4人の例(荘島 2012 の p.7, p.11): 1 = 好き，7 = 嫌い --------------------------
ex3 <- matrix(c(NA, 1, 7,
                7, NA, 1,
                1, 7, NA), 3, byrow = TRUE, dimnames = list(LETTERS[1:3], LETTERS[1:3]))
ex4a <- matrix(c(NA, 1, 1, 7,
                 1, NA, 1, 7,
                 7, 7, NA, 1,
                 1, 1, 7, NA), 4, byrow = TRUE, dimnames = list(LETTERS[1:4], LETTERS[1:4]))
ex4b <- matrix(c(NA, 7, 7, 7,
                 1, NA, 1, 7,
                 1, 7, NA, 1,
                 1, 1, 7, NA), 4, byrow = TRUE, dimnames = list(LETTERS[1:4], LETTERS[1:4]))

# 4. 千野 (1997, p.13) のソシオメトリックデータ(荘島 2012 の改変版。小さいほど好意) -----------------
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

examples <- list(ex3 = ex3, ex4a = ex4a, ex4b = ex4b, chino = chino)
fits0 <- map(examples, vonmises_mds_try, lambda = 0)
fits1 <- map(examples, vonmises_mds_try, lambda = 1)
tibble(data = names(examples),
       stress_l0 = map_dbl(fits0, "stress"), stress_l1 = map_dbl(fits1, "stress"),
       conv_l0 = map_lgl(fits0, "converged"), conv_l1 = map_lgl(fits1, "converged"))

out_png <- file.path(Sys.getenv("VM_OUT", tempdir()), "vonmises_try.png")
png(out_png, width = 1600, height = 800, res = 110)
par(mfrow = c(2, 4), mar = c(4, 4, 2, 1))
iwalk(fits0, ~ plot_vm_try(.x, main = paste(.y, "lambda = 0")))
iwalk(fits1, ~ plot_vm_try(.x, main = paste(.y, "lambda = 1")))
dev.off()
out_png
