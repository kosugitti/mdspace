# 楕円モデル(ellipse_mds)のテスト: 書籍表4.4 との照合・勾配・円モデルとの関係 ---------------

occ <- as.matrix(utils::read.csv(test_path("fixtures", "okada_imaizumi_2024_table1_19.csv"),
                                 row.names = 1, check.names = FALSE))
use <- row(occ) != col(occ) & !is.na(occ)
cells_of <- function(model, alpha = -1, beta = 1) {
  list(s = occ[use], j = row(occ)[use], k = col(occ)[use], n = 10, decreasing = TRUE,
       model = model, alpha = alpha, beta = beta, kappa = Inf, axis_max = Inf)
}

# 書籍表4.4(p.118)の2次元解: 座標2列と半軸2列(行の順はデータと同じ)
tab44 <- matrix(c(1.191, 0.656, 1.665, 0.624, 1.025, 0.047, 1.290, 0.394,
                  0.988, -1.020, 1.165, 1.154, -0.188, 0.395, 1.143, 0.599,
                  -0.311, -0.342, 0.149, 0.175, 0.215, -0.689, 0.262, 0.117,
                  -0.453, 0.956, 1.175, 0.821, -1.241, 0.355, 1.072, 1.303,
                  -0.703, -0.697, 0.679, 0.430, -0.523, 0.339, 1.091, 0.628),
                10, byrow = TRUE)

test_that("S of the book's Table 4.4 solution is 0.324", {
  par <- mdspace:::.ell_par(tab44[, 1:2], tab44[, 3:4], "II")
  expect_equal(mdspace:::.ell_stress(par, cells_of("II"), 2), 0.324, tolerance = 5e-4 / 0.324)
})

test_that("the analytic gradient agrees with finite differences", {
  num_grad <- function(par, cells, p) {
    h <- 1e-6
    vapply(seq_along(par), function(i) {
      e <- replace(numeric(length(par)), i, h)
      (mdspace:::.ell_stress(par + e, cells, p) - mdspace:::.ell_stress(par - e, cells, p)) / (2 * h)
    }, numeric(1))
  }
  set.seed(1)
  for (ab in list(c(-1, 1), c(1, -1), c(-1, -1), c(1, 1), c(0, 1), c(-1, 0))) {
    cl <- cells_of("II", ab[1], ab[2])
    par <- c(rnorm(20), log(runif(20, 0.3, 1.5)))
    expect_equal(mdspace:::.ell_grad(par, cl, 2), num_grad(par, cl, 2), tolerance = 1e-5)
    cl <- cells_of("III", ab[1], ab[2])
    par <- c(rnorm(20), runif(10, 0, 0.5))
    expect_equal(mdspace:::.ell_grad(par, cl, 2), num_grad(par, cl, 2), tolerance = 1e-5)
  }
})

test_that("Model IIIB reproduces circle_mds()", {
  fe <- ellipse_mds(occ, ndim = 2, model = "III", nrandom = 3)
  fc <- circle_mds(occ, ndim = 2, nrandom = 3, kappa = Inf)
  expect_equal(fe$stress, fc$stress, tolerance = 1e-4)
  expect_s3_class(fe, "mdspace_ellipse")
  expect_length(fe$radius, 10)
  expect_equal(min(fe$radius), 0)
})

test_that("Model II minimizes S below the book and warns about degenerate ellipses", {
  expect_warning(
    expect_warning(fit <- ellipse_mds(occ, ndim = 2, nrandom = 2, check_orientation = FALSE),
                   "long and thin"),
    "very large")
  expect_lt(fit$stress, 0.324)
  # 正規化: 重心が原点・座標の二乗和が n・次元は分散の大きい順
  expect_equal(unname(colMeans(fit$conf)), c(0, 0), tolerance = 1e-8)
  expect_equal(sum(fit$conf^2), 10, tolerance = 1e-8)
  expect_true(all(diff(colSums(fit$conf^2)) <= 0))
  expect_true(all(fit$semiaxis > 0))
  # 返した布置と半軸から計算し直した S が一致する
  par <- mdspace:::.ell_par(fit$conf, fit$semiaxis, "II")
  expect_equal(mdspace:::.ell_stress(par, cells_of("II"), 2), fit$stress, tolerance = 1e-6)
})

test_that("bounds kappa and axis_max are respected and reported", {
  fit <- suppressWarnings(ellipse_mds(occ, ndim = 2, nrandom = 1, kappa = 3, axis_max = 1.2,
                                      check_orientation = FALSE))
  A <- fit$semiaxis
  expect_true(all(apply(A, 1, max) / apply(A, 1, min) <= 3 + 1e-6))
  expect_true(all(exp(rowMeans(log(A))) <= 1.2 + 1e-6))
  expect_true(all(fit$at_bound))
  msgs <- character(0)
  withCallingHandlers(
    ellipse_mds(occ, ndim = 2, nrandom = 0, kappa = 3, axis_max = 1.2, check_orientation = FALSE),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("aspect ratio.*reached the bound", msgs)))
  expect_true(any(grepl("size.*reached the bound", msgs)))
})

test_that("the orientation check returns S for each angle and direction", {
  fit <- suppressWarnings(ellipse_mds(occ, ndim = 2, nrandom = 0, kappa = 3, axis_max = 1.2,
                                      angles = c(10, 20)))
  o <- fit$orientation
  expect_equal(nrow(o), 2)
  expect_named(o, c("dims", "angle", "counterclockwise", "clockwise"))
  # 半軸を回して座標だけ合わせ直した S は，回さない解より小さくならない(わずかな数値誤差は許す)
  expect_true(all(c(o$counterclockwise, o$clockwise) >= fit$stress - 1e-3))
})

test_that("missing cells are allowed", {
  x <- occ
  x[1, 2] <- x[5, 7] <- x[9, 3] <- NA
  fit <- suppressWarnings(ellipse_mds(x, ndim = 2, nrandom = 1, check_orientation = FALSE))
  expect_true(is.na(fit$mjk[1, 2]))
  expect_false(is.na(fit$mjk[2, 1]))
  expect_true(is.finite(fit$stress))
})

test_that("in one dimension Model II and Model III fit the same model", {
  # 1次元では楕円の半軸は1つだけ＝円と同じ m_jk になる
  par2 <- mdspace:::.ell_par(matrix(tab44[, 1]), matrix(tab44[, 3]), "II")
  par3 <- mdspace:::.ell_par(matrix(tab44[, 1]), tab44[, 3], "III")
  expect_equal(mdspace:::.ell_parts(par2, cells_of("II"), 1)$m,
               mdspace:::.ell_parts(par3, cells_of("III"), 1)$m)
})

test_that("invalid arguments are rejected", {
  expect_error(ellipse_mds(occ, model = "I"), "not implemented")
  expect_error(ellipse_mds(occ, alpha = 0, beta = 0), "not both 0")
  expect_error(ellipse_mds(occ, alpha = 2), "-1, 0 or 1")
  expect_error(ellipse_mds(occ, kappa = 0.5), "at least 1")
  expect_error(ellipse_mds(occ, axis_max = 0), "positive")
})

test_that("print, summary and plot methods work", {
  fit <- suppressWarnings(ellipse_mds(occ, ndim = 2, nrandom = 0, kappa = 3, axis_max = 1.2,
                                      check_orientation = FALSE))
  expect_output(print(fit))
  expect_no_error(summary(fit))
  pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(plot(fit))
  expect_no_error(plot(fit, which = "semiaxis"))
})
