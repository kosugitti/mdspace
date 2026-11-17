# 二相三元の円モデル(circle_mds にリストを渡す)のテスト: 表1.14と nASYMS の出力 ---------------

raw <- utils::read.csv(test_path("fixtures", "okada_imaizumi_2024_table1_14.csv"),
                       check.names = FALSE)
nations <- lapply(1:2, function(i) {
  m <- as.matrix(raw[raw$source == i, -(1:2)])
  rownames(m) <- colnames(m)
  m
})
ns_conf <- utils::read.csv(test_path("fixtures", "nasyms_twomode_nations_conf.csv"))
ns_wt <- utils::read.csv(test_path("fixtures", "nasyms_twomode_nations_weights.csv"))
cells_of <- function(Xs) {
  lapply(Xs, function(X) {
    use <- row(X) != col(X) & !is.na(X)
    list(s = X[use], j = row(X)[use], k = col(X)[use], decreasing = TRUE)
  })
}

test_that("S and w reproduce nASYMS for its own solutions (2 to 5 dimensions)", {
  cells <- cells_of(nations)
  for (q in 2:5) {
    d <- ns_conf[ns_conf$ndim == q, ]
    ww <- ns_wt[ns_wt$ndim == q, ]
    X <- as.matrix(d[, 2 + seq_len(q)])
    U <- as.matrix(ww[, paste0("u", seq_len(q))])
    S <- mdspace:::.c3_stress(c(X, d$radius, U), 10, q, 2, cells)
    expect_equal(S, ww$S[1], tolerance = 2e-3 / ww$S[1])
    w <- vapply(1:2, function(i) {
      o <- mdspace:::.c3_msrc(X, d$radius, U[i, ], cells[[i]])
      sqrt(10 / sum((o$m - mean(o$m))^2))
    }, numeric(1))
    expect_equal(w, ww$w, tolerance = 2e-3)
  }
})

test_that("the analytic gradient agrees with finite differences", {
  cells <- cells_of(nations)
  set.seed(1)
  par <- c(rnorm(20), runif(10, 0, 0.5), runif(4, 0.5, 1.5))
  h <- 1e-6
  num <- vapply(seq_along(par), function(a) {
    e <- replace(numeric(length(par)), a, h)
    (mdspace:::.c3_stress(par + e, 10, 2, 2, cells) -
       mdspace:::.c3_stress(par - e, 10, 2, 2, cells)) / (2 * h)
  }, numeric(1))
  expect_equal(mdspace:::.c3_grad(par, 10, 2, 2, cells), num, tolerance = 1e-5)
})

test_that("minimized S is no worse than nASYMS", {
  fit <- suppressWarnings(circle_mds(nations, ndim = 2, ndim_max = 5))
  expect_s3_class(fit, "mdspace_circle")
  expect_true(isTRUE(fit$threeway))
  ref <- c("5" = 0.400, "4" = 0.431, "3" = 0.473, "2" = 0.540)
  expect_true(all(fit$path <= ref + 5e-4))
  expect_lte(fit$stress, 0.540)
  expect_equal(fit$stress, sqrt(mean(fit$stress_source^2)), tolerance = 1e-8)
  expect_equal(names(fit$weights), c("source", "w", "u1", "u2"))
  expect_true(all(fit$aspect <= 10 + 1e-8))
  expect_equal(sum(fit$conf^2), 10, tolerance = 1e-8)
  expect_equal(min(fit$radius), 0)
})

test_that("a true ellipse model is recovered", {
  set.seed(5)
  n <- 15
  X <- matrix(rnorm(n * 2), n, 2)
  r <- runif(n, 0, 0.6)
  U <- rbind(c(1.5, 0.75), c(0.8, 1.6), c(1, 1))
  Xs <- lapply(1:3, function(i) {
    dif_j <- row(diag(n)); dif_k <- col(diag(n))
    M <- matrix(NA_real_, n, n)
    for (j in 1:n) for (k in 1:n) if (j != k) {
      dx <- X[j, ] - X[k, ]
      d <- sqrt(sum(dx^2)); v <- d / sqrt(sum((dx / U[i, ])^2))
      M[j, k] <- d - v * (r[j] - r[k])
    }
    -M
  })
  fit <- circle_mds(Xs, ndim = 2, nrandom = 3)
  expect_lt(fit$stress, 1e-4)
  expect_equal(unname(fit$aspect), c(2, 2, 1), tolerance = 0.02)
  expect_gt(cor(fit$radius, r), 0.99)
  # 長軸の向き(ソース1は第1次元，ソース2は第2次元)は回転で入れ替わりうるので，両者が逆なことを見る
  long <- apply(as.matrix(fit$weights[, c("u1", "u2")]), 1, which.max)
  expect_true(long[1] != long[2])
})

test_that("bad input is an error and print works", {
  expect_error(circle_mds(list(nations[[1]])), "at least two")
  expect_error(circle_mds(list(nations[[1]], nations[[2]][1:5, 1:5])), "same size")
  fit <- suppressWarnings(circle_mds(nations, ndim = 2, nrandom = 1))
  expect_output(print(fit), "Two-mode three-way")
  pdf(NULL)
  expect_invisible(plot(fit))
  dev.off()
})
