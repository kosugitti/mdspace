# HFM と非対称分解のテスト ------------------------------------------------------

# jscp_lesson2.Rmd(2019)の 4x4 の例
x4 <- matrix(c(0, 1, 1, 7, 1, 0, 1, 7, 7, 7, 0, 1, 1, 1, 7, 0), 4, byrow = TRUE,
             dimnames = list(LETTERS[1:4], LETTERS[1:4]))

test_that("asym_decomp splits X into S + K", {
  dec <- asym_decomp(x4)
  expect_equal(unname(dec$S + dec$K), unname(x4))
  expect_true(isSymmetric(unname(dec$S)))
  expect_equal(unname(dec$K), -t(unname(dec$K)))
  expect_equal(unname(dec$ss["total"]), sum(x4^2))
})

test_that("hfm_mds reproduces H with all dimensions", {
  fit <- hfm_mds(x4, ndim = 4)
  Z <- fit$coord
  Hhat <- Z %*% diag(sign(fit$eigen)) %*% Conj(t(Z))
  expect_equal(Re(Hhat), unname(asym_decomp(x4)$S), ignore_attr = TRUE)
  expect_equal(Im(Hhat), unname(asym_decomp(x4)$K), ignore_attr = TRUE)
  expect_equal(fit$stress, 0, tolerance = 1e-12)
  expect_equal(sum(fit$gof), 1)
})

test_that("hfm_mds eigenvalues agree with the 2019 lesson code", {
  # 2019 年の自作関数と同じ手順(エルミート化→eigen→絶対値の順)
  H <- matrix(complex(real = (x4 + t(x4)) / 2, imaginary = (x4 - t(x4)) / 2), 4)
  ev <- eigen(H)$values
  ev <- ev[order(abs(ev), decreasing = TRUE)]
  expect_equal(hfm_mds(x4, ndim = 4)$eigen, Re(ev))
})

test_that("std rotation puts the principal axis on the real axis", {
  fit <- hfm_mds(x4, ndim = 2, rotate = "std")
  for (t in 1:2) {
    z <- fit$coord[, t]
    expect_equal(Im(sum(z^2)), 0, tolerance = 1e-10)
    expect_gte(Re(sum(z^2)), 0)
    # どの角度に回しても Re の二乗和はこれより大きくならない
    best <- max(sapply(seq(0, pi, length.out = 721),
                       function(a) sum(Re(z * exp(1i * a))^2)))
    expect_lte(best, sum(Re(z)^2) + 1e-8)
  }
})

test_that("varimax rotation maximizes the variance of Re(z^2)", {
  fit <- hfm_mds(x4, ndim = 2, rotate = "varimax")
  crit <- function(z) stats::var(Re(z^2))
  for (t in 1:2) {
    z <- fit$coord[, t]
    best <- max(sapply(seq(0, pi, length.out = 721),
                       function(a) crit(z * exp(1i * a))))
    expect_lte(best, crit(z) + 1e-8)
  }
})

test_that("rotation does not change the reproduced H", {
  a <- hfm_mds(x4, ndim = 4, rotate = "varimax")
  b <- hfm_mds(x4, ndim = 4, rotate = "none")
  ha <- a$coord %*% diag(sign(a$eigen)) %*% Conj(t(a$coord))
  hb <- b$coord %*% diag(sign(b$eigen)) %*% Conj(t(b$coord))
  expect_equal(ha, hb, ignore_attr = TRUE)
})

test_that("a symmetric matrix gives no imaginary parts", {
  s <- (x4 + t(x4)) / 2
  fit <- hfm_mds(s, ndim = 2)
  expect_equal(Im(fit$coord), matrix(0, 4, 2), ignore_attr = TRUE, tolerance = 1e-12)
})

test_that("the returned object has the common mdspace slots", {
  fit <- hfm_mds(x4, ndim = 2)
  expect_s3_class(fit, c("mdspace_hfm", "mdspace"))
  expect_equal(dim(fit$conf), c(4, 4))
  expect_equal(fit$labels, LETTERS[1:4])
  expect_error(hfm_mds(x4[1:3, ]), "square")
  expect_error(hfm_mds(x4, ndim = 5), "ndim")
})
