# 欠測の型の集計と MINERS 法のテスト ------------------------------------------------

# 階数 2 のエルミート行列から作った，欠測を正確に戻せるはずの行列
make_lowrank <- function(n = 6, seed = 1) {
  set.seed(seed)
  U <- qr.Q(qr(matrix(complex(real = rnorm(n * 2), imaginary = rnorm(n * 2)), n, 2)))
  H <- U %*% diag(c(20, 8)) %*% Conj(t(U))
  X <- Re(H) + Im(H)
  dimnames(X) <- list(LETTERS[1:n], LETTERS[1:n])
  return(X)
}

test_that("missing_pattern classifies the four types", {
  X <- make_lowrank(6)
  X[1, 2] <- NA                      # 片側
  X[3, 4] <- NA; X[4, 3] <- NA       # 両側
  X[5, 5] <- NA                      # 対角
  p <- missing_pattern(X)
  expect_equal(p$one_side, 1)
  expect_equal(p$both_sides, 2)
  expect_equal(p$diagonal, 1)
  expect_equal(p$absent, 0)
  Y <- make_lowrank(6)
  Y[2, -2] <- NA; Y[-2, 2] <- NA     # いない人(対角は観測されていても absent)
  p2 <- missing_pattern(list(a = X, b = Y))
  expect_equal(p2$absent, c(0, 1))
  expect_equal(p2$absent_objects[2], "B")
  expect_equal(p2$both_sides[2], 0)
})

test_that("one-sided cells are filled by the transposed value", {
  X <- make_lowrank(6)
  X0 <- X
  X[1, 2] <- NA
  fit <- hfm_mds(X, ndim = 2, missing = "impute")
  expect_equal(fit$data[1, 2], X0[2, 1])
  expect_equal(fit$imputed$type, "one_side")
})

test_that("MINERS recovers diagonal and two-sided cells of a low-rank matrix", {
  X <- make_lowrank(8, seed = 2)
  X0 <- X
  X[3, 3] <- NA
  X[2, 6] <- NA; X[6, 2] <- NA
  fit <- hfm_mds(X, ndim = 2, missing = "impute", tol = 1e-10, maxit = 5000)
  expect_equal(fit$data[3, 3], X0[3, 3], tolerance = 1e-4)
  expect_equal(fit$data[2, 6], X0[2, 6], tolerance = 1e-4)
  expect_equal(fit$data[6, 2], X0[6, 2], tolerance = 1e-4)
  expect_setequal(fit$imputed$type, c("diagonal", "both_sides"))
  expect_gt(fit$impute_niter, 1)
})

test_that("missing values are an error unless imputed; absent objects are always an error", {
  X <- make_lowrank(6)
  X[1, 2] <- NA
  expect_error(hfm_mds(X), "impute")
  Y <- make_lowrank(6)
  Y[2, -2] <- NA; Y[-2, 2] <- NA
  expect_error(hfm_mds(Y, missing = "impute"), "absent")
  expect_error(gstatis_mds(list(Y, make_lowrank(6)), missing = "impute"), "absent")
})

test_that("gstatis_mds imputes each source and reports the cells", {
  X1 <- make_lowrank(6, 3)
  X2 <- make_lowrank(6, 4)
  X2[1, 1] <- NA
  X2[2, 3] <- NA
  fit <- gstatis_mds(list(T1 = X1, T2 = X2), ndim = 2, missing = "impute")
  expect_equal(nrow(fit$imputed), 2)
  expect_equal(unique(fit$imputed$source), "T2")
  expect_false(anyNA(fit$data$T2))
})
