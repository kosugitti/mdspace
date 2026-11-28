# 等高線マップモデル contour_mds() のテスト -------------------------------------

# 小杉・藤原 (2004) の調査対象者 C と同じ評定・広島県の布置(研究：等高線/ContourMapping/test/CL広島S4.csv)
make_baseball <- function() {
  conf <- matrix(c(2.0407, -0.4754, -0.2996, 0.4323, 1.0778, 1.0394,
                   -1.6636, 0.5439, -0.4076, -0.5291, -0.7478, -1.0111),
                 ncol = 2, byrow = TRUE,
                 dimnames = list(c("Giants", "Dragons", "Tigers", "Carp", "Swallows", "Baystars"),
                                 NULL))
  list(conf = conf, raw = c(Giants = 1, Dragons = 2, Tigers = 1, Carp = 5, Swallows = 2,
                            Baystars = 2))
}

# 千野版 Abelson.map()(develop/Abelson.map.R)をそのまま写したもの
abelson_map_chino <- function(dat, locations) {
  X <- dat[, 1]
  Y <- dat[, 2]
  P <- dat[, 3]
  un <- matrix(1, nrow(locations), 1, )
  Xs <- un %*% X
  Ys <- un %*% Y
  dm <- ((locations[, 1] - Xs)^2 + (locations[, 2] - Ys)^2) + 1
  V <- t(P %*% (1 / t(dm)))
  return(as.vector(V))
}

test_that("heights at the objects follow equation (1)", {
  b <- make_baseball()
  fit <- contour_mds(b$conf, b$raw, center = 3)
  expect_s3_class(fit, c("mdspace_contour", "mdspace"))
  expect_equal(unname(fit$height[1, ]),
               c(-2.737, -2.011, -2.845, 0.794, -2.361, -2.031), tolerance = 1e-3)
  expect_equal(fit$discrepancy, fit$ratings - fit$height)
  # ずれ = ほかの対象から来る力の和の符号を変えたもの
  W <- 1 / (1 + as.matrix(dist(b$conf))^2)
  diag(W) <- 0
  expect_equal(unname(fit$discrepancy[1, ]), -unname(as.vector(W %*% (b$raw - 3))))
})

test_that("valences on a grid agree with Chino's Abelson.map()", {
  b <- make_baseball()
  fit <- contour_mds(b$conf, b$raw, center = 3)
  loc <- as.matrix(expand.grid(seq(-3, 3, by = 0.5), seq(-2, 2, by = 0.5)))
  expect_equal(unname(predict(fit, loc)[, 1]),
               abelson_map_chino(cbind(b$conf, b$raw - 3), loc))
  expect_equal(predict(fit), t(fit$height))
})

test_that("valences are invariant to rotation and translation and linear in ratings", {
  b <- make_baseball()
  th <- 0.7
  Q <- matrix(c(cos(th), sin(th), -sin(th), cos(th)), 2)
  conf2 <- b$conf %*% Q + matrix(c(3, -1), 6, 2, byrow = TRUE)
  rownames(conf2) <- rownames(b$conf)
  f1 <- contour_mds(b$conf, b$raw, center = 3)
  f2 <- contour_mds(conf2, b$raw, center = 3)
  expect_equal(f1$height, f2$height)
  f3 <- contour_mds(b$conf, 2 * (b$raw - 3))
  expect_equal(f3$height, 2 * f1$height)
})

test_that("ratings are matched by name and several persons are allowed", {
  b <- make_baseball()
  r <- rbind(A = b$raw, B = rev(b$raw))
  colnames(r) <- names(b$raw)
  r2 <- r[, rev(colnames(r))]
  f1 <- contour_mds(b$conf, r, center = 3)
  f2 <- contour_mds(b$conf, r2, center = 3)
  expect_equal(f1$height, f2$height)
  expect_equal(dim(f1$height), c(2L, 6L))
  expect_equal(f1$height["A", ], contour_mds(b$conf, b$raw, center = 3)$height[1, ])
  bad <- b$raw
  names(bad)[1] <- "Kyojin"
  expect_error(contour_mds(b$conf, bad), "do not match")
  expect_error(contour_mds(b$conf, unname(b$raw)[-1]), "6 values")
})

test_that("normalize scales the configuration to unit mean squared distance", {
  b <- make_baseball()
  fit <- contour_mds(b$conf * 5, b$raw, center = 3, normalize = TRUE)
  expect_equal(sum(fit$conf^2) / 6, 1)
  expect_equal(unname(colMeans(fit$conf)), c(0, 0))
  expect_equal(fit$height, contour_mds(b$conf, b$raw, center = 3, normalize = TRUE)$height)
})

test_that("missing ratings carry no charge and three dimensions work", {
  b <- make_baseball()
  r <- b$raw
  r["Carp"] <- NA
  fit <- contour_mds(b$conf, r, center = 3)
  r0 <- b$raw - 3
  r0["Carp"] <- 0
  expect_equal(unname(fit$height[1, ]),
               unname(as.vector((1 / (1 + as.matrix(dist(b$conf))^2)) %*% r0)))
  expect_true(is.na(fit$discrepancy[1, "Carp"]))
  expect_error(contour_mds(b$conf, rep(NA_real_, 6)), "at least one")
  conf3 <- cbind(b$conf, c(0.5, -0.2, 0.1, 0.3, -0.4, -0.3))
  f3 <- contour_mds(conf3, b$raw, center = 3)
  expect_equal(unname(f3$height[1, ]),
               unname(as.vector((1 / (1 + as.matrix(dist(conf3))^2)) %*% (b$raw - 3))))
  expect_error(contour_mds(b$conf[, 1, drop = FALSE], b$raw), "at least 2")
})

test_that("methods print, summarize and plot", {
  b <- make_baseball()
  fit <- contour_mds(b$conf, b$raw, center = 3)
  expect_output(print(fit), "contour")
  expect_output(print(summary(fit)), "discrepancy")
  fit2 <- contour_mds(b$conf, rbind(A = b$raw, B = rev(b$raw)), center = 3)
  expect_output(print(summary(fit2)), "Discrepancy")
  pdf(NULL)
  on.exit(dev.off())
  g <- plot(fit)
  expect_equal(dim(g$z), c(101L, 101L))
  expect_silent(plot(fit2, person = "B", type = "persp"))
  expect_error(plot(fit2, person = "C"), "person")
})
