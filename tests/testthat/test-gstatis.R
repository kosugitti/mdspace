# GSTATIS のテスト ---------------------------------------------------------------

x4 <- matrix(c(0, 1, 1, 7, 1, 0, 1, 7, 7, 7, 0, 1, 1, 1, 7, 0), 4, byrow = TRUE,
             dimnames = list(LETTERS[1:4], LETTERS[1:4]))

test_that("identical occasions give RV = 1 and the HFM of each occasion", {
  fit <- gstatis_mds(list(x4, x4, x4), ndim = 2)
  expect_equal(unname(fit$rv), matrix(1, 3, 3))
  # 重みは GRV の固有ベクトルそのもの(1/sqrt(3) が3つ)
  expect_equal(unname(fit$weights[, 1]), rep(1 / sqrt(3), 3))
  h <- hfm_mds(x4, ndim = 2)
  expect_equal(fit$eigen, h$eigen * sqrt(3))
  # 時点がすべて同じなら軌跡は1点に重なり，Lp ノルムは 0
  for (t in 2:3) {
    expect_equal(fit$trajectory[, , t], fit$trajectory[, , 1])
  }
  expect_equal(max(fit$lpnorm), 0)
  expect_length(fit$occasion, 3)
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
  expect_gt(sum(fit$weights[, 1]), 0)
  expect_equal(sum(fit$weights[, 1]^2), 1)
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

test_that("trajectories follow the revised step 3 (D U* h)", {
  set.seed(4)
  xs <- lapply(1:3, function(t) x4 + matrix(runif(16), 4))
  fit <- gstatis_mds(xs, ndim = 3, project = "row")
  cp <- fit$compromise[[1]]
  U <- cp$vectors[, 1:3]
  D <- ifelse(cp$eigen[1:3] < 0, -1i, 1)
  H2 <- .hermitian(xs[[2]])
  # 対象 1 の行ベクトルの像
  img <- as.vector(Conj(t(U)) %*% H2[1, ]) * D
  expect_equal(unname(fit$trajectory[1, , 2]), img)
  # Lp ノルムは像の差の長さ
  expect_equal(unname(fit$lpnorm[1, 2, 1]),
               sqrt(sum(Mod(fit$trajectory[1, , 1] - fit$trajectory[1, , 2])^2)))
  # 列で射影すると t(H) の行になる
  fitc <- gstatis_mds(xs, ndim = 3, project = "column")
  imgc <- as.vector(Conj(t(U)) %*% H2[, 1]) * D
  expect_equal(unname(fitc$trajectory[1, , 2]), imgc)
})
