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

test_that("impute_ndim must be smaller than n", {
  X <- make_lowrank(6)
  X[3, 3] <- NA
  expect_error(hfm_mds(X, ndim = 6, missing = "impute"), "impute_ndim")
  expect_error(hfm_mds(X, ndim = 2, impute_ndim = 6, missing = "impute"), "impute_ndim")
})

test_that("initial values are the mean of nonzero observed cells", {
  Y <- matrix(c(NA, 5, 0, 0, 1,
                4, NA, 0, 0, 2,
                0, 0, NA, 6, 3,
                0, 1, 7, NA, 0,
                2, 1, 2, 0, NA), 5, byrow = TRUE,
              dimnames = list(letters[1:5], letters[1:5]))
  # maxit = 0 なら反復せず初期値のまま返る: 初期値を直接確かめる
  imp <- suppressWarnings(mdspace:::.impute_miners(Y, p = 1, maxit = 0))
  expect_equal(unname(diag(imp$data)), c(3, 2.6, 4.5, 14 / 3, 11 / 6))
  # 行と列がすべてゼロなら観測値の平均(= 0)に戻る
  Z <- Y
  Z[1, -1] <- 0; Z[-1, 1] <- 0
  imp2 <- suppressWarnings(mdspace:::.impute_miners(Z, p = 1, maxit = 0))
  expect_equal(unname(imp2$data[1, 1]), 0)
})

test_that("stopping at the first iteration warns", {
  X <- make_lowrank(6)
  X[3, 3] <- NA
  # 階数 5(= n - 1) なら再構成がほぼ今の行列のまま
  expect_warning(mdspace:::.impute_miners(X, p = 5, tol = 1e-2), "first iteration")
})

# 円・楕円・フォンミーゼスのモデルの欠測(セルを当てはめから外す) ---------------------------

# 真の楕円モデル(二相三元)から作った，誤差のない類似度のソース3つ
make_twomode <- function(n = 15, seed = 5) {
  set.seed(seed)
  X <- matrix(rnorm(n * 2), n, 2)
  r <- runif(n, 0, 0.6)
  U <- rbind(c(1.5, 0.75), c(0.8, 1.6), c(1, 1))
  Xs <- lapply(1:3, function(i) {
    M <- matrix(NA_real_, n, n)
    for (j in 1:n) for (k in 1:n) if (j != k) {
      dx <- X[j, ] - X[k, ]
      d <- sqrt(sum(dx^2))
      v <- d / sqrt(sum((dx / U[i, ])^2))
      M[j, k] <- d - v * (r[j] - r[k])
    }
    dimnames(M) <- list(LETTERS[1:n], LETTERS[1:n])
    -M
  })
  names(Xs) <- c("T1", "T2", "T3")
  list(Xs = Xs, r = r)
}

test_that("two-mode circle model: missing cells are left out and the solution is still recovered", {
  g <- make_twomode()
  Xs <- g$Xs
  Xs$T1["A", "B"] <- NA                       # 片側
  Xs$T2["C", "D"] <- Xs$T2["D", "C"] <- NA    # 両側
  Xs$T3["E", "I"] <- NA
  fit <- circle_mds(Xs, ndim = 2, nrandom = 3)
  expect_lt(fit$stress, 1e-4)
  expect_gt(cor(fit$radius, g$r), 0.99)
  expect_equal(unname(fit$aspect), c(2, 2, 1), tolerance = 0.02)
  expect_false(anyNA(fit$conf))
})

test_that("two-mode circle model: an object absent at one source is placed from the others", {
  g <- make_twomode()
  Xs <- g$Xs
  Xs$T2["G", ] <- NA
  Xs$T2[, "G"] <- NA
  expect_equal(missing_pattern(Xs)$absent_objects, c("", "G", ""))
  fit <- circle_mds(Xs, ndim = 2, nrandom = 3)
  expect_lt(fit$stress, 1e-4)
  expect_gt(cor(fit$radius, g$r), 0.99)
  expect_true(all(is.finite(fit$stress_source)))
  # 誤差のないデータなので，欠けのない解と点間距離が一致する
  full <- circle_mds(g$Xs, ndim = 2, nrandom = 3)
  expect_gt(cor(c(dist(fit$conf)), c(dist(full$conf))), 0.999)
})

test_that("two-mode circle model: an object absent at every source, or an empty source, is an error", {
  g <- make_twomode()
  Xs <- lapply(g$Xs, function(m) {
    m["G", ] <- NA
    m[, "G"] <- NA
    m
  })
  expect_error(circle_mds(Xs, ndim = 2, nrandom = 0), "G have no observed proximities in any matrix")
  Ys <- g$Xs
  Ys$T3[] <- NA
  expect_error(circle_mds(Ys, ndim = 2, nrandom = 0), "Source\\(s\\) T3 have no observed")
})

test_that("one-mode models stop on an object with no observed proximities", {
  set.seed(1)
  n <- 10
  z <- matrix(rnorm(n * 2), n, 2)
  D <- as.matrix(dist(z)) + matrix(runif(n * n, 0, 0.2), n, n)
  dimnames(D) <- list(LETTERS[1:n], LETTERS[1:n])
  D["D", ] <- NA
  D[, "D"] <- NA
  expect_error(circle_mds(D, ndim = 2, type = "dissimilarity"), "D have no observed proximities")
  expect_error(ellipse_mds(D, ndim = 2, type = "dissimilarity"), "D have no observed proximities")
  expect_error(vonmises_mds(D, type = "dissimilarity"), "D have no observed proximities")
  # 片側だけでも観測があれば止めない(行か列のどちらかで決まる)
  E <- as.matrix(dist(z)) + matrix(runif(n * n, 0, 0.2), n, n)
  E[4, ] <- NA
  expect_no_error(suppressWarnings(circle_mds(E, ndim = 2, type = "dissimilarity")))
})
