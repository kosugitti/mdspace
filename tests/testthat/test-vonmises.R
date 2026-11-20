# フォンミーゼスモデル(vonmises_mds)のテスト: 勾配・識別の枠・パラメタリカバリ・返り値 ----------------

# 千野の10対象の非類似度(develop/20261118_vonmises_map_check.R と同じ)
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

# 内部の cells を作る(vonmises_mds() の中と同じ手順)
vm_cells <- function(G, free = TRUE, prior = list()) {
  pr <- list(x_sd = 1, kappa_shape = 2, kappa_rate = 0.5, c_a = 2, c_b = 4)
  pr[names(prior)] <- prior
  diag(G) <- NA
  G <- G / sqrt(mean(G^2, na.rm = TRUE))
  use <- !is.na(G)
  list(n = nrow(G), i = row(G)[use], j = col(G)[use], g = G[use], free = free,
       kmax = if (free) 20 else 6.5, prior = pr)
}

# モデルから非類似度を作る: xi_ij = (1 - pi_ij) delta_ij，pi_ij = c_i exp{kappa_i (cos(theta_ij - mu_i) - 1)}
vm_simulate <- function(X, mu, kappa, c, sd = 0) {
  n <- nrow(X)
  G <- matrix(NA_real_, n, n)
  for (i in seq_len(n)) for (j in seq_len(n)) if (i != j) {
    d <- X[j, ] - X[i, ]
    p <- c[i] * exp(kappa[i] * (cos(atan2(d[2], d[1]) - mu[i]) - 1))
    G[i, j] <- (1 - p) * sqrt(sum(d^2))
  }
  G + sd * matrix(stats::rnorm(n * n), n)
}

test_that("the analytic gradient agrees with finite differences", {
  num_grad <- function(par, cells) {
    h <- 1e-6
    vapply(seq_along(par), function(k) {
      e <- replace(numeric(length(par)), k, h)
      (mdspace:::.vm_nlp(par + e, cells) - mdspace:::.vm_nlp(par - e, cells)) / (2 * h)
    }, numeric(1))
  }
  G <- chino
  G[2, 5] <- G[7, 1] <- NA
  set.seed(1)
  for (free in c(TRUE, FALSE)) {
    cells <- vm_cells(G, free)
    par <- mdspace:::.vm_random(cells)
    expect_equal(mdspace:::.vm_gr(par, cells), num_grad(par, cells), tolerance = 1e-5)
  }
})

test_that("the log posterior does not depend on rotation, reflection or translation", {
  cells <- vm_cells(chino)
  set.seed(2)
  X <- matrix(rnorm(20), 10, 2)
  mu <- runif(10, -pi, pi)
  kappa <- runif(10, 0.5, 4)
  cc <- runif(10, 0.1, 0.6)
  base <- mdspace:::.vm_nlp(mdspace:::.vm_pack(X, mu, kappa, cc, 0.5, cells), cells)
  # 回転(角度 a)と平行移動: 向き mu も同じだけ回る
  a <- 0.7
  R <- matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2)
  Xr <- sweep(X %*% R, 2, c(3, -2), "+")
  expect_equal(mdspace:::.vm_nlp(mdspace:::.vm_pack(Xr, mu - a, kappa, cc, 0.5, cells), cells), base)
  # 鏡映(第2軸を反転): 向きも反転する
  Xm <- cbind(X[, 1], -X[, 2])
  expect_equal(mdspace:::.vm_nlp(mdspace:::.vm_pack(Xm, -mu, kappa, cc, 0.5, cells), cells), base)
  # pack -> unpack で点間距離が保たれる
  p <- mdspace:::.vm_unpack(mdspace:::.vm_pack(X, mu, kappa, cc, 0.5, cells), cells)
  expect_equal(as.matrix(dist(p$X)), as.matrix(dist(X)), ignore_attr = TRUE)
  expect_equal(p$kappa, kappa)
  expect_equal(p$c, cc)
})

test_that("parameters are recovered from noiseless data", {
  set.seed(3)
  n <- 8
  X <- matrix(rnorm(2 * n), n, 2)
  mu <- runif(n, -pi, pi)
  kappa <- runif(n, 1, 4)
  cc <- runif(n, 0.2, 0.6)
  G <- vm_simulate(X, mu, kappa, cc)
  fit <- vonmises_mds(G, nrandom = 5)
  # 布置は大きさだけ自由(非類似度の尺度で割っている)＝点間距離の相関で見る
  expect_gt(stats::cor(c(dist(fit$conf)), c(dist(X))), 0.99)
  expect_lt(fit$stress, 0.05)
  # 向きは回転・鏡映をそろえてから比べる
  A <- sweep(fit$conf, 2, colMeans(fit$conf))
  B <- sweep(X, 2, colMeans(X))
  s <- svd(t(B) %*% A)
  R <- s$u %*% t(s$v)
  v <- cbind(cos(mu), sin(mu)) %*% R
  err <- abs(((atan2(v[, 2], v[, 1]) - fit$mu + pi) %% (2 * pi)) - pi)
  expect_lt(max(err), 10 * pi / 180)
})

test_that("the returned solution is normalized and internally consistent", {
  for (om in c("free", "fixed")) {
    fit <- vonmises_mds(chino, omega = om, nrandom = 3)
    expect_s3_class(fit, "mdspace_vonmises")
    # 重心が原点・主軸(交差積が対角)
    expect_equal(unname(colMeans(fit$conf)), c(0, 0), tolerance = 1e-8)
    cp <- crossprod(fit$conf)
    expect_lt(abs(cp[1, 2]), 1e-8)
    expect_gte(cp[1, 1], cp[2, 2])
    expect_true(all(fit$mu >= -pi & fit$mu < pi))
    expect_true(all(fit$kappa > 0 & fit$kappa < if (om == "free") 20 else 6.5))
    expect_true(all(fit$c > 0 & fit$c < 1))
    if (om == "fixed") {
      expect_equal(unname(fit$c), 1 / (2 * pi * besselI(unname(fit$kappa), 0, expon.scaled = TRUE)))
    }
    # 最良の初期値の対数事後密度を返す
    expect_equal(fit$logpost, max(fit$starts))
    expect_equal(fit$logpost, fit$starts[fit$best_start])
    # 返した布置・向き・kappa・c から xi と pi を計算し直すと一致する(xi は元の単位)
    d <- outer(seq_len(10), seq_len(10), Vectorize(function(i, j) {
      v <- fit$conf[j, ] - fit$conf[i, ]
      p <- fit$c[i] * exp(fit$kappa[i] * (cos(atan2(v[2], v[1]) - fit$mu[i]) - 1))
      (1 - p) * sqrt(sum(v^2))
    }))
    diag(d) <- NA
    expect_equal(fit$xi, d * fit$scale, ignore_attr = TRUE, tolerance = 1e-8)
    g <- chino / fit$scale
    expect_equal(fit$stress, sqrt(sum((g - d)^2, na.rm = TRUE) / sum(g^2, na.rm = TRUE)),
                 tolerance = 1e-8)
    expect_length(fit$anchor, 3)
  }
})

test_that("similarities are reversed into dissimilarities", {
  sim <- max(chino, na.rm = TRUE) + min(chino, na.rm = TRUE) - chino
  f1 <- vonmises_mds(chino, nrandom = 1)
  f2 <- vonmises_mds(sim, type = "similarity", nrandom = 1)
  expect_equal(f2$logpost, f1$logpost)
  expect_equal(f2$conf, f1$conf)
})

test_that("missing cells are allowed", {
  G <- chino
  G[1, 2] <- G[5, 7] <- G[9, 3] <- NA
  fit <- vonmises_mds(G, nrandom = 1)
  expect_true(is.na(fit$xi[1, 2]))
  expect_false(is.na(fit$xi[2, 1]))
  expect_true(is.finite(fit$logpost))
})

test_that("approximate standard errors are returned with hessian = TRUE", {
  fit <- vonmises_mds(chino, omega = "fixed", nrandom = 1, hessian = TRUE)
  expect_named(fit$se, c("mu", "kappa", "c"))
  expect_equal(nrow(fit$se), 10)
  expect_true(all(is.finite(as.matrix(fit$se)) & as.matrix(fit$se) >= 0))
})

test_that("the random number state of the caller is kept", {
  set.seed(10)
  a <- runif(1)
  set.seed(10)
  invisible(vonmises_mds(chino, nrandom = 2, seed = 5))
  expect_equal(runif(1), a)
})

test_that("invalid input is rejected", {
  expect_error(vonmises_mds(matrix(c(NA, 1, 1, NA), 2)), "At least 3")
  G <- chino
  G[1, 2] <- -1
  expect_error(vonmises_mds(G), "nonnegative")
})

test_that("print, summary and plot methods work", {
  fit <- vonmises_mds(chino, nrandom = 1, hessian = TRUE)
  expect_output(print(fit), "omega")
  s <- summary(fit)
  expect_true(all(c("mu_deg", "kappa", "c", "se_mu_deg") %in% colnames(s$conf)))
  pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(plot(fit))
})
