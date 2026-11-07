# GSTATIS のテスト ---------------------------------------------------------------

x4 <- matrix(c(0, 1, 1, 7, 1, 0, 1, 7, 7, 7, 0, 1, 1, 1, 7, 0), 4, byrow = TRUE,
             dimnames = list(LETTERS[1:4], LETTERS[1:4]))

test_that("identical occasions give RV = 1 and the HFM of each occasion", {
  fit <- gstatis_mds(list(x4, x4, x4), ndim = 2)
  expect_equal(unname(fit$rv), matrix(1, 3, 3))
  expect_equal(unname(fit$weights[, 1]), rep(1 / 3, 3))
  h <- hfm_mds(x4, ndim = 2)
  expect_equal(fit$eigen, h$eigen)
  expect_equal(fit$coord, h$coord, ignore_attr = TRUE)
  # 各時点の軌跡は妥協の座標と一致する
  for (t in 1:3) {
    expect_equal(fit$trajectory[, , t], fit$coord, ignore_attr = TRUE)
  }
})

test_that("RV coefficient matches the 2001 definition", {
  set.seed(2)
  xs <- lapply(1:3, function(t) x4 + matrix(runif(16), 4))
  fit <- gstatis_mds(xs, ndim = 2)
  H <- lapply(xs, function(m) {
    matrix(complex(real = (m + t(m)) / 2, imaginary = (m - t(m)) / 2), 4)
  })
  # VB の coefficientHFM: sum_k sum_l A(k,l,i) * A(l,k,j)
  tr <- function(a, b) Re(sum(diag(a %*% b)))
  rv12 <- tr(H[[1]], H[[2]]) / sqrt(tr(H[[1]], H[[1]]) * tr(H[[2]], H[[2]]))
  expect_equal(unname(fit$rv[1, 2]), rv12)
  expect_true(isSymmetric(unname(fit$rv)))
  expect_equal(sum(fit$weights[, 1]), 1)
})

test_that("array input, ncomp, threshold and scale work", {
  set.seed(3)
  arr <- array(sample(1:7, 4 * 4 * 4, TRUE), c(4, 4, 4))
  fit <- gstatis_mds(arr, ndim = 1, ncomp = 2)
  expect_s3_class(fit, c("mdspace_gstatis", "mdspace"))
  expect_length(fit$compromise, 2)
  expect_equal(dim(fit$trajectory), c(4, 1, 4))
  fit2 <- gstatis_mds(arr, ndim = 1, ncomp = NULL, threshold = 0)
  expect_length(fit2$compromise, 4)
  fit3 <- gstatis_mds(arr, ndim = 1, scale = "norm")
  expect_equal(fit3$rv, fit$rv)
  expect_error(gstatis_mds(list(x4)), "at least two")
  expect_error(gstatis_mds(list(x4, x4[1:3, 1:3])), "same size")
})
