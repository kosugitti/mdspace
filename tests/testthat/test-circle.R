# 円モデル(circle_mds)のテスト: 書籍表1.19と nASYMS の解との照合 ---------------------

occ <- as.matrix(utils::read.csv(test_path("fixtures", "okada_imaizumi_2024_table1_19.csv"),
                                 row.names = 1, check.names = FALSE))
nasyms_all <- utils::read.csv(test_path("fixtures", "nasyms_radius_occupation_all.csv"))
# nASYMS が表示した S(2026-11-17 の実行)
nasyms_S <- c(0.499, 0.372, 0.276, 0.202, 0.183)
nasyms <- nasyms_all[nasyms_all$ndim == 2, ]
names(nasyms)[3:4] <- c("dim1", "dim2")

test_that("S reproduces nASYMS for its own solutions in 1 to 5 dimensions", {
  use <- row(occ) != col(occ) & !is.na(occ)
  cells <- list(s = occ[use], j = row(occ)[use], k = col(occ)[use], n = 10,
                decreasing = TRUE)
  for (p in 1:5) {
    d <- nasyms_all[nasyms_all$ndim == p, ]
    Z <- as.matrix(d[, 2 + seq_len(p)])
    # 布置は小数3桁なので S も 3桁目の丸めの範囲で一致
    expect_equal(mdspace:::.circle_stress(c(Z, d$radius), cells, p), nasyms_S[p],
                 tolerance = 1e-3 / nasyms_S[p])
  }
})

test_that("the analytic gradient agrees with finite differences", {
  set.seed(1)
  use <- row(occ) != col(occ) & !is.na(occ)
  cells <- list(s = occ[use], j = row(occ)[use], k = col(occ)[use], n = 10,
                decreasing = TRUE)
  par <- c(rnorm(20), runif(10))
  h <- 1e-6
  num <- vapply(seq_along(par), function(i) {
    e <- replace(numeric(length(par)), i, h)
    (mdspace:::.circle_stress(par + e, cells, 2) -
       mdspace:::.circle_stress(par - e, cells, 2)) / (2 * h)
  }, numeric(1))
  expect_equal(mdspace:::.circle_grad(par, cells, 2), num, tolerance = 1e-5)
})

test_that("minimized S is no worse than nASYMS from 5 down to 1 dimensions", {
  fit <- circle_mds(occ, ndim = 1, ndim_max = 5, nrandom = 0)
    ref <- stats::setNames(rev(nasyms_S), 5:1)
  expect_true(all(fit$path <= ref + 5e-4))
  expect_s3_class(fit, "mdspace_circle")
})

test_that("the 2-dimensional solution matches nASYMS up to rotation", {
  fit <- circle_mds(occ, ndim = 2, ndim_max = 5, nrandom = 0)
  expect_equal(fit$stress, 0.372, tolerance = 2e-3)
  # 正規化: 重心 0，座標の二乗和 n，最小の半径 0
  expect_equal(unname(colMeans(fit$conf)), c(0, 0), tolerance = 1e-8)
  expect_equal(sum(fit$conf^2), 10, tolerance = 1e-8)
  expect_equal(min(fit$radius), 0)
  # 直交回転(と鏡映)をそろえて座標と半径を比べる
  Z <- as.matrix(nasyms[, c("dim1", "dim2")])
  sv <- svd(t(fit$conf) %*% Z)
  rot <- fit$conf %*% sv$u %*% t(sv$v)
  expect_lt(max(abs(rot - Z)), 0.1)
  expect_lt(max(abs(fit$radius - nasyms$radius)), 0.1)
  # m と disparity の対応
  expect_equal(sum(is.na(fit$mjk)), 10)
  expect_equal(fit$stress,
               sqrt(sum((fit$mjk - fit$disparity)^2, na.rm = TRUE) /
                      sum((fit$mjk - mean(fit$mjk, na.rm = TRUE))^2, na.rm = TRUE)),
               tolerance = 1e-8)
})

test_that("dissimilarities give the same solution as negated similarities", {
  a <- circle_mds(occ, ndim = 2, nrandom = 0)
  b <- circle_mds(-occ, ndim = 2, type = "dissimilarity", nrandom = 0)
  expect_equal(a$stress, b$stress, tolerance = 1e-6)
})

test_that("random starts are kept only when better, and methods run", {
  fit <- circle_mds(occ, ndim = 2, nrandom = 2, seed = 1)
  expect_length(fit$starts, 3)
  expect_equal(fit$stress, min(fit$starts))
  expect_output(print(fit), "circle")
  expect_output(print(summary(fit)), "radius")
  pdf(NULL)
  expect_invisible(plot(fit))
  dev.off()
  expect_error(circle_mds(occ, ndim = 10), "ndim")
})

test_that("random starts are on by default, reproducible, and keep the RNG state", {
  set.seed(99)
  before <- .Random.seed
  a <- circle_mds(occ, ndim = 2)
  expect_identical(.Random.seed, before)
  expect_length(a$starts, 11)
  b <- circle_mds(occ, ndim = 2)
  expect_equal(a$stress, b$stress)
  expect_lte(a$stress, a$starts[1])
})

test_that("dimensions beyond the positive eigenvalues still move", {
  # 誤差のない 20 対象 10 次元: Torgerson の正の固有値が 10 個に足りない例
  set.seed(1)
  n <- 20; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  r <- runif(n)
  S <- -(as.matrix(dist(X)) - outer(r, r, "-"))
  diag(S) <- NA
  fit <- expect_no_warning(circle_mds(S, ndim = p, nrandom = 0))
  expect_true(all(apply(fit$conf, 2, sd) > 0.01))
  # 乱数の初期値を足せば誤差のない真の解(S = 0)に着く
  expect_lt(circle_mds(S, ndim = p)$stress, 1e-6)
})

test_that("off-diagonal missing cells are left out of the fit", {
  # 誤差のない 10 対象 2 次元に欠測を入れても真の解(S = 0)に着く
  set.seed(3)
  n <- 10
  X <- matrix(rnorm(n * 2), n, 2)
  r <- runif(n, 0, 0.5)
  S <- -(as.matrix(dist(X)) - outer(r, r, "-"))
  diag(S) <- NA
  S[1, 2] <- NA          # 片側
  S[3, 4] <- S[4, 3] <- NA  # 両側
  fit <- circle_mds(S, ndim = 2)
  expect_lt(fit$stress, 1e-6)
  expect_true(is.na(fit$mjk[1, 2]) && is.na(fit$disparity[3, 4]))
  expect_false(is.na(fit$mjk[2, 1]))
})
