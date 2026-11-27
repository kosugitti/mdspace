# S3 メソッド(print・summary・plot)がクラスごとにそろっているかのテスト -----------------

x4 <- matrix(c(0, 1, 1, 7, 1, 0, 1, 7, 7, 7, 0, 1, 1, 1, 7, 0), 4, byrow = TRUE,
             dimnames = list(LETTERS[1:4], LETTERS[1:4]))
raw <- utils::read.csv(test_path("fixtures", "okada_imaizumi_2024_table1_14.csv"),
                       check.names = FALSE)
nations <- lapply(1:2, function(i) {
  m <- as.matrix(raw[raw$source == i, -(1:2)])
  rownames(m) <- colnames(m)
  m
})
names(nations) <- c("first", "second")

# 図を捨てる装置を開いて式を評価する
with_null_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  force(code)
}

test_that("every mdspace class has print, summary and plot methods", {
  classes <- c("mdspace", "mdspace_hfm", "mdspace_gstatis", "mdspace_circle",
               "mdspace_ellipse", "mdspace_vonmises", "mdspace_asym")
  ns <- asNamespace("mdspace")
  for (cl in classes) {
    for (generic in c("print", "summary", "plot")) {
      # 子クラスの print は共通の print.mdspace を継承してもよいが，ここでは全クラスが自前で持つ
      expect_true(exists(paste0(generic, ".", cl), envir = ns, inherits = FALSE),
                  info = paste0(generic, ".", cl))
    }
  }
})

test_that("the common methods work on a bare mdspace object", {
  conf <- matrix(c(1, -1, 0, 0, 0, 0, 1, -1), 4, dimnames = list(LETTERS[1:4], NULL))
  obj <- mdspace:::.new_mdspace("demo", conf = conf, stress = 0.1234, niter = 5L)
  expect_s3_class(obj, c("mdspace_demo", "mdspace"), exact = TRUE)
  expect_output(print(obj), "Stress: 0.1234")
  utils::capture.output(expect_invisible(print(obj)))
  s <- summary(obj)
  expect_s3_class(s, "summary.mdspace")
  expect_output(print(s), "Configuration")
  expect_invisible(with_null_device(plot(obj)))
  # 1 次元の布置も描ける
  one <- mdspace:::.new_mdspace("demo", conf = conf[, 1, drop = FALSE])
  expect_invisible(with_null_device(plot(one)))
  # extra の表は名前つきで順に出る(データフレームは数値の列だけ丸める)
  s$extra <- list("A table" = data.frame(name = "a", value = 1.23456),
                  "A vector" = c(v = 2.34567))
  expect_output(print(s), "A table:")
  expect_output(print(s), "1.2346")
  expect_output(print(s), "2.3457")
})

test_that("hfm methods: print, summary, plot", {
  fit <- hfm_mds(x4, ndim = 2)
  expect_output(print(fit), "Unexplained SS")
  utils::capture.output(expect_invisible(print(fit)))
  s <- summary(fit)
  expect_s3_class(s, "summary.mdspace")
  expect_true(is.na(s$stress))
  tab <- s$extra[["Eigenvalues and goodness of fit"]]
  expect_equal(unname(tab["eigenvalue", ]), fit$eigen[1:2])
  expect_equal(unname(tab["cumulative", 2]), 1 - fit$unexplained)
  expect_output(print(s), "Eigenvalues and goodness of fit")
  expect_output(print(s), "Unexplained SS")
  expect_invisible(with_null_device(plot(fit)))
  expect_invisible(with_null_device(plot(fit, dim = 2)))
})

test_that("gstatis methods: print, summary, plot", {
  set.seed(2)
  xs <- lapply(1:3, function(t) x4 + matrix(runif(16), 4))
  names(xs) <- c("t1", "t2", "t3")
  fit <- gstatis_mds(xs, ndim = 2, ncomp = 2)
  expect_output(print(fit), "RV eigenvalues")
  utils::capture.output(expect_invisible(print(fit)))
  s <- summary(fit)
  expect_s3_class(s, "summary.mdspace")
  w <- s$extra[["Weights of the sources (GRV eigenvectors)"]]
  expect_equal(dim(w), c(3L, 2L))
  expect_equal(rownames(w), c("t1", "t2", "t3"))
  expect_equal(unname(w), unname(fit$weights))
  d <- s$extra[["Distance of each source from the compromise"]]
  expect_equal(unname(d[, 1]), unname(fit$compromise[[1]]$distance))
  expect_output(print(s), "RV eigenvalues")
  expect_output(print(s), "First compromise")
  for (type in c("trajectory", "compromise", "interstructure")) {
    expect_invisible(with_null_device(plot(fit, type = type)))
  }
  expect_invisible(with_null_device(plot(fit, objects = c("A", "C"), scale = "max")))
})

test_that("asym_decomp methods: print, summary, plot", {
  dec <- asym_decomp(x4)
  expect_output(print(dec), "Share of asymmetry")
  utils::capture.output(expect_invisible(print(dec)))
  s <- summary(dec)
  expect_s3_class(s, "summary.mdspace_asym")
  expect_equal(unname(s$ss["share", "symmetric"] + s$ss["share", "skew"]), 1)
  expect_equal(s$asymmetry, dec$asymmetry)
  # 歪対称部分の行和は対象全体で 0，絶対値の行和は |K| の総和
  expect_equal(sum(s$objects[, "skew_sum"]), 0)
  expect_equal(sum(s$objects[, "abs_skew_sum"]), sum(abs(dec$K)))
  expect_equal(rownames(s$objects), LETTERS[1:4])
  expect_output(print(s), "Row sums of the skew-symmetric part")
  utils::capture.output(expect_invisible(print(s)))
  expect_invisible(with_null_device(plot(dec)))
})

test_that("summary of the two-mode three-way circle model shows the sources", {
  fit <- suppressWarnings(circle_mds(nations, ndim = 2, nrandom = 1, seed = 1))
  s <- summary(fit)
  expect_equal(colnames(s$conf), c("D1", "D2", "radius"))
  wt <- s$extra[[1]]
  expect_equal(wt$source, c("first", "second"))
  expect_equal(names(wt), c("source", "w", "u1", "u2", "aspect", "S"))
  expect_equal(wt$S, unname(fit$stress_source))
  expect_output(print(s), "Sources")
  # 一相二元の summary に extra は付かない
  one <- circle_mds(nations[[1]], ndim = 2, nrandom = 0)
  expect_null(summary(one)$extra)
})

test_that("plot of the two-mode three-way circle model draws the ellipses of sources", {
  fit <- suppressWarnings(circle_mds(nations, ndim = 2, nrandom = 1, seed = 1))
  # 描いた線の数と軸の範囲を記録する(graphics::lines を数える代わりに記録装置を使う)
  draw <- function(...) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    grDevices::dev.control(displaylist = "enable")
    plot(fit, ...)
    list(usr = graphics::par("usr"), rec = grDevices::recordPlot())
  }
  common <- draw()
  s1 <- draw(source = 1)
  s2 <- draw(source = "second")
  all <- draw(source = "all")
  # 軸の範囲はどの描き方でも同じ(並べて比べられる)
  expect_equal(s1$usr, common$usr)
  expect_equal(s2$usr, common$usr)
  expect_equal(all$usr, common$usr)
  # 全ソースの楕円が枠に収まる
  U <- as.matrix(fit$weights[, c("u1", "u2")])
  far <- max(abs(fit$conf) + outer(fit$radius, apply(U, 2, max)))
  expect_true(far <= max(abs(common$usr)))
  expect_invisible(with_null_device(plot(fit, source = 2, labels = FALSE)))
  # ソースの指定の誤り
  expect_error(with_null_device(plot(fit, source = 3)), "does not match")
  expect_error(with_null_device(plot(fit, source = "third")), "does not match")
  expect_error(with_null_device(plot(fit, source = 1:2)), "single")
  # 一相二元の解に source は渡せない
  one <- circle_mds(nations[[1]], ndim = 2, nrandom = 0)
  expect_error(with_null_device(plot(one, source = 1)), "two-mode three-way")
  expect_invisible(with_null_device(plot(one)))
})

test_that("plot of the two-mode three-way circle model works in 1 and 3 dimensions", {
  fit1 <- suppressWarnings(circle_mds(nations, ndim = 1, nrandom = 0))
  expect_invisible(with_null_device(plot(fit1, source = "all")))
  fit3 <- suppressWarnings(circle_mds(nations, ndim = 3, nrandom = 0))
  expect_invisible(with_null_device(plot(fit3, dims = c(1, 3), source = 1)))
})

test_that("plot of the asymmetry weights of the two-mode three-way circle model", {
  fit <- suppressWarnings(circle_mds(nations, ndim = 2, nrandom = 1, seed = 1))
  usr <- with_null_device({
    plot(fit, which = "weights")
    graphics::par("usr")
  })
  # 全ソースの重みが枠に収まり，原点が入る
  U <- as.matrix(fit$weights[, c("u1", "u2")])
  expect_true(max(U) <= min(usr[c(2, 4)]))
  expect_true(usr[1] <= 0 && usr[3] <= 0)
  expect_invisible(with_null_device(plot(fit, which = "weights", labels = FALSE)))
  fit3 <- suppressWarnings(circle_mds(nations, ndim = 3, nrandom = 0))
  expect_invisible(with_null_device(plot(fit3, which = "weights", dims = c(1, 3))))
  # 1 次元と一相二元では描けない
  fit1 <- suppressWarnings(circle_mds(nations, ndim = 1, nrandom = 0))
  expect_error(with_null_device(plot(fit1, which = "weights")), "two or more dimensions")
  one <- circle_mds(nations[[1]], ndim = 2, nrandom = 0)
  expect_error(with_null_device(plot(one, which = "weights")), "two-mode three-way")
})
