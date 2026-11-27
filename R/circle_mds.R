# 円モデル(radius-distance model): Okada & Imaizumi (1987) -----------------------

#' Nonmetric asymmetric MDS by the radius-distance (circle) model
#'
#' Each object `j` is a point `x_j` in a `p`-dimensional Euclidean space with a
#' circle (sphere, hypersphere) of radius `r_j` centered at it. The
#' asymmetric proximity from `j` to `k` is represented by
#' `m_jk = d_jk - r_j + r_k`, the distance from the near edge of the circle of
#' `j` to the far edge of the circle of `k`, which is assumed to be
#' monotonically decreasing in similarities (increasing in dissimilarities).
#' The fit is measured by
#' `S = sqrt(sum((m_jk - mhat_jk)^2) / sum((m_jk - mbar)^2))`
#' (Kruskal's stress formula 2 applied to `m_jk`), where `mhat_jk` are the
#' disparities obtained by monotone regression (primary approach to ties).
#' `S` is minimized over the coordinates and radii.
#'
#' The initial configuration is "rational" as in Okada and Imaizumi (1987):
#' Torgerson's method applied to the symmetrized dissimilarities, with the
#' additive constant that makes them Euclidean (Cailliez, 1983, as in
#' [stats::cmdscale()] with `add = TRUE`; Okada & Imaizumi, 1997, add
#' Torgerson's constant), gives the configuration in `ndim_max` dimensions, and the least squares solution of
#' `r_k - r_j = a_jk` (`a_jk` the skew-symmetric part of the dissimilarities)
#' gives the radii. The solution in `ndim_max` dimensions is then reduced one
#' dimension at a time by its principal axes down to `ndim`. Because the
#' radius-distance model has local minima, `nrandom` random starts (10 by
#' default) are also tried in `ndim` dimensions and the solution with the
#' smallest `S` is returned (Okada & Imaizumi, 2024, Section 4.1.1).
#'
#' With few objects the solution may not be unique: a simulation with
#' `n = 5` (20 off-diagonal cells) reached `S = 0` even with error while the
#' coordinates differed from the true ones by up to 0.8 (normalized scale), and
#' one case with `n = 10` in three dimensions did the same. The configuration
#' should then be interpreted with care. See
#' `develop/20261117_circle_recovery.R` in the source repository.
#'
#' The solution is normalized so that the configuration is centered with the
#' sum of squared coordinates equal to `n`, and the smallest radius is 0 (the
#' radii are determined only up to an additive constant), as in the nASYMS
#' program.
#'
#' **Two-mode three-way model.** When `x` is a list of `N` matrices (sources),
#' the model of Okada and Imaizumi (1997) is fitted: a common configuration
#' of points and circles, a symmetry weight `w_i` and asymmetry weights `u_it`
#' for each source. In source `i` the circle of object `j` becomes an ellipse
#' with semiaxes `u_it r_j`, and
#' `m_jki = w_i d_jk - v_jki (r_j - r_k)` with
#' `v_jki = w_i d_jk / sqrt(sum_t ((x_jt - x_kt) / u_it)^2)`.
#' `S = sqrt(mean_i S_i^2)`, where `S_i` is stress formula 2 within source
#' `i`. The asymmetry weights are fixed to have geometric mean 1 over all
#' sources and dimensions (they and the radii are otherwise determined only
#' up to reciprocal scaling), the aspect ratio of each source is bounded by
#' `kappa`, and the problem is solved by L-BFGS-B. `w_i` does not affect `S`
#' and is set afterwards so that `sum_jk (m_jki - mbar_i)^2 = n`, as in the
#' original paper and nASYMS. The initial configuration is the Torgerson
#' solution (with the additive constant) of the mean of the symmetrized
#' dissimilarities; when reducing from `ndim_max`, the dimensions with the
#' smallest contribution are dropped. The result additionally has
#' `weights` (data frame with `w` and `u1`, `u2`, ...), `aspect`,
#' `stress_source` and `at_bound`.
#'
#' The minimization uses the analytic gradient of `S` with the disparities
#' held fixed, which is the gradient of `S` as a function of the coordinates
#' and radii (Kruskal, 1964b), and [stats::optim()] (BFGS) instead of the
#' steepest descent of the original program.
#'
#' @param x A square matrix of asymmetric proximities, or a list of such
#'   matrices (one per source) for the two-mode three-way model. The diagonal
#'   is ignored; other missing cells are left out of the fit.
#' @param ndim Number of dimensions.
#' @param type `"similarity"` (default) or `"dissimilarity"`. If `x` has a
#'   `type` attribute set by [as_proximity()], that is used.
#' @param ndim_max Dimensionality at which the rational initial configuration
#'   is built; the solutions are reduced from there down to `ndim`.
#' @param nrandom Number of random starts tried in addition to the rational
#'   start (`0` uses the rational start only).
#' @param maxit,reltol Passed to [stats::optim()].
#' @param kappa Two-mode three-way model only: the upper bound of the ratio
#'   of the longest to the shortest semiaxis of the ellipses of a source
#'   (`max_t u_it / min_t u_it`). Without it the ellipses can collapse to line
#'   segments while the radii diverge. A warning is given when a solution
#'   reaches the bound.
#' @param seed Random seed for the random starts. The state of the random
#'   number generator is restored afterwards. `NULL` uses (and advances) the
#'   current state of the generator.
#' @return An object of class `c("mdspace_circle", "mdspace")` with `conf`
#'   (n x ndim), `radius`, `stress` (`S`), `niter` (gradient evaluations of
#'   the chosen solution), `mjk` and `disparity` (n x n matrices, `NA` on the
#'   diagonal and missing cells), `path` (`S` at each dimensionality from
#'   `ndim_max` to `ndim` for the rational start), `starts` (`S` of each
#'   start, the rational one first), `best_start` (which start was chosen),
#'   and `converged`.
#' @references Okada, A., & Imaizumi, T. (1987). Nonmetric multidimensional
#'   scaling of asymmetric proximities. *Behaviormetrika*, 14(21), 81-96.
#'
#' Okada, A., & Imaizumi, T. (2024). *Applied multidimensional scaling of
#'   asymmetric relationships*. Springer.
#'
#' Cailliez, F. (1983). The analytical solution of the additive constant
#'   problem. *Psychometrika*, 48(2), 305-308.
#'
#' Okada, A., & Imaizumi, T. (1997). Asymmetric multidimensional scaling of
#'   two-mode three-way proximities. *Journal of Classification*, 14, 195-224.
#'
#' Kruskal, J. B. (1964b). Nonmetric multidimensional scaling: A numerical
#'   method. *Psychometrika*, 29(2), 115-129.
#' @examples
#' x <- matrix(c(NA, 5, 2, 1,
#'               3, NA, 4, 2,
#'               1, 2, NA, 6,
#'               2, 1, 4, NA), 4, byrow = TRUE)
#' fit <- circle_mds(x, ndim = 1)
#' fit
#' @export
circle_mds <- function(x, ndim = 2, type = c("similarity", "dissimilarity"),
                       ndim_max = ndim, nrandom = 10, maxit = 1000,
                       reltol = 1e-10, seed = 1, kappa = 10) {
  cl <- match.call()
  # 行列のリストなら二相三元(Okada & Imaizumi, 1997) ------------------------------------
  if (is.list(x) && !is.data.frame(x)) {
    type <- match.arg(type)
    return(.circle3(x, ndim, type, ndim_max, nrandom, maxit, seed, kappa, cl))
  }
  if (!is.null(attr(x, "type")) && missing(type)) {
    type <- attr(x, "type")
  }
  type <- match.arg(type)
  X <- as_proximity(x, type = type)
  n <- nrow(X)
  if (ndim < 1 || ndim >= n) {
    stop("`ndim` must be between 1 and n - 1.")
  }
  if (ndim_max < ndim || ndim_max >= n) {
    stop("`ndim_max` must be between `ndim` and n - 1.")
  }

  # 当てはめに使うセル(対角と欠測を除く) ---------------------------------------
  use <- row(X) != col(X) & !is.na(X)
  cells <- list(s = X[use], j = row(X)[use], k = col(X)[use], n = n,
                decreasing = type == "similarity")
  if (length(cells$s) < 3) {
    stop("Too few observed off-diagonal cells.")
  }

  # 合理的初期布置: 最大次元で作り，主軸で1次元ずつ下ろす -------------------------
  init <- .circle_init(X, ndim_max, type)
  Z <- init$conf
  r <- init$radius
  path <- numeric(0)
  fit <- NULL
  for (p in seq(ndim_max, ndim)) {
    if (ncol(Z) > p) {
      Z <- .principal_axes(Z, p)
    }
    fit <- .circle_optim(Z, r, cells, maxit, reltol)
    Z <- fit$conf
    r <- fit$radius
    path[as.character(p)] <- fit$stress
  }
  best <- fit
  best_start <- 1L
  starts <- fit$stress

  # 乱数の初期値(局所解を避けるため既定で10本) -----------------------------------------
  if (nrandom > 0) {
    if (!is.null(seed)) {
      # 呼び出し側の乱数の状態を壊さないよう，終わったら元に戻す
      old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
      on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()), add = TRUE)
      set.seed(seed)
    }
    for (i in seq_len(nrandom)) {
      Z0 <- matrix(stats::rnorm(n * ndim), n, ndim)
      r0 <- stats::runif(n, 0, 0.5)
      f <- .circle_optim(Z0, r0, cells, maxit, reltol)
      starts <- c(starts, f$stress)
      if (f$stress < best$stress) {
        best <- f
        best_start <- i + 1L
      }
    }
  }

  # 正規化して返す ---------------------------------------------------------------
  nr <- .circle_normalize(best$conf, best$radius)
  conf <- nr$conf
  rownames(conf) <- rownames(X)
  colnames(conf) <- paste0("D", seq_len(ndim))
  radius <- stats::setNames(nr$radius, rownames(X))
  m <- .circle_m(conf, radius)
  M <- D <- matrix(NA_real_, n, n, dimnames = dimnames(X))
  M[use] <- m[use]
  D[use] <- .monotone(M[use], cells$s, cells$decreasing)

  out <- .new_mdspace("circle", conf = conf, stress = best$stress,
                      niter = best$niter, data = X, call = cl, ndim = ndim,
                      radius = radius, mjk = M, disparity = D, type = type,
                      path = path, starts = starts, best_start = best_start, converged = best$converged)
  return(out)
}

# 初期布置と初期半径(内部) ---------------------------------------------------------
# 類似度は最大値から引いて非類似度にし，対称部に Torgerson 法(stats::cmdscale)，
# 歪対称部 a_jk = r_k - r_j の最小二乗解(和 0)を半径の初期値にする
.circle_init <- function(X, p, type, additive = TRUE) {
  # 加える定数で Torgerson 法の結果が変わらないよう，非対角の最小を 0 にそろえる
  delta <- if (type == "similarity") -X else X
  diag(delta) <- NA
  delta <- delta - min(delta, na.rm = TRUE)
  diag(delta) <- 0
  # 欠測は転置の値，それもなければ全体の平均で埋めて初期値だけに使う
  miss <- is.na(delta)
  delta[miss] <- t(delta)[miss]
  delta[is.na(delta)] <- mean(delta, na.rm = TRUE)
  sym <- (delta + t(delta)) / 2
  skew <- (delta - t(delta)) / 2
  # 加算定数(Okada & Imaizumi 1997 の初期値は Torgerson 1952 の加算定数を足す)。
  # ここでは stats::cmdscale(add = TRUE) の Cailliez (1983) の定数＝非対角に足すと
  # ユークリッド距離になる最小の定数を使う
  if (additive) {
    conf <- suppressWarnings(stats::cmdscale(sym, k = p, add = TRUE)$points)
  } else {
    conf <- suppressWarnings(stats::cmdscale(sym, k = p))
  }
  conf <- matrix(conf, nrow = nrow(X))
  # 正の固有値が p 個に足りないときは，残りの次元を固有ベクトルの向きに小さく与える
  # (0 のままだとその次元の勾配が常に 0 になり，動かなくなる)
  if (ncol(conf) < p) {
    n <- nrow(X)
    J <- diag(n) - 1 / n
    ev <- eigen(-0.5 * J %*% sym^2 %*% J, symmetric = TRUE)
    # 中心化で必ず出る定数ベクトル(布置に差を与えない)は除く
    cand <- setdiff(which(abs(colSums(ev$vectors)) < 1e-6), seq_len(ncol(conf)))
    extra <- cand[seq_len(p - ncol(conf))]
    small <- 0.01 * sqrt(max(ev$values[1], 1e-12))
    conf <- cbind(conf, ev$vectors[, extra, drop = FALSE] * small)
  }
  radius <- colMeans(skew)
  return(list(conf = conf, radius = radius - min(radius)))
}

# 布置を主軸の上位 p 本に射影する(内部) ---------------------------------------------
.principal_axes <- function(Z, p) {
  Z <- sweep(Z, 2, colMeans(Z))
  V <- svd(Z)$v[, seq_len(p), drop = FALSE]
  return(Z %*% V)
}

# 正規化: 重心を原点に，座標の二乗和を n に，最小の半径を 0 に(内部) ----------------------
.circle_normalize <- function(Z, r) {
  Z <- sweep(Z, 2, colMeans(Z))
  a <- sqrt(nrow(Z) / sum(Z^2))
  r <- r * a
  return(list(conf = Z * a, radius = r - min(r)))
}

# m_jk = d_jk - r_j + r_k の行列(内部) -------------------------------------------------
.circle_m <- function(Z, r) {
  D <- as.matrix(stats::dist(Z))
  return(D - outer(r, r, "-"))
}

# 単調回帰(第1の方法: 同順位の中は m の順に並べ替えてよい)(内部) ------------------------
.monotone <- function(m, s, decreasing) {
  o <- if (decreasing) order(-s, m) else order(s, m)
  mh <- numeric(length(m))
  mh[o] <- stats::isoreg(m[o])$yf
  return(mh)
}

# S とその勾配(内部) -------------------------------------------------------------------
.circle_parts <- function(par, cells, nd) {
  p <- nd
  n <- cells$n
  Z <- matrix(par[seq_len(n * p)], n, p)
  r <- par[n * p + seq_len(n)]
  dif <- Z[cells$j, , drop = FALSE] - Z[cells$k, , drop = FALSE]
  d <- sqrt(rowSums(dif^2))
  m <- d - r[cells$j] + r[cells$k]
  return(list(dif = dif, d = d, m = m))
}

.circle_stress <- function(par, cells, nd) {
  q <- .circle_parts(par, cells, nd)
  mh <- .monotone(q$m, cells$s, cells$decreasing)
  return(sqrt(sum((q$m - mh)^2) / sum((q$m - mean(q$m))^2)))
}

.circle_grad <- function(par, cells, nd) {
  n <- cells$n
  p <- nd
  q <- .circle_parts(par, cells, p)
  m <- q$m
  mh <- .monotone(m, cells$s, cells$decreasing)
  A <- sum((m - mh)^2)
  B <- sum((m - mean(m))^2)
  S <- sqrt(A / B)
  if (S < 1e-12) {
    return(numeric(length(par)))
  }
  # dS/dm_jk
  g <- ((m - mh) - S^2 * (m - mean(m))) / (S * B)
  # 重なった点(d = 0)では距離の勾配を 0 とする
  w <- ifelse(q$d > 1e-12, g / q$d, 0)
  idx <- factor(c(cells$j, cells$k), levels = seq_len(n))
  GZ <- matrix(0, n, p)
  for (t in seq_len(p)) {
    v <- w * q$dif[, t]
    GZ[, t] <- as.vector(tapply(c(v, -v), idx, sum, default = 0))
  }
  gr <- as.vector(tapply(c(-g, g), idx, sum, default = 0))
  return(c(GZ, gr))
}

# 1つの初期値からの最小化(内部) -----------------------------------------------------
# (optim の引数 par に部分一致しないよう，次元数は nd という名前で渡す)
.circle_optim <- function(Z, r, cells, maxit, reltol) {
  p <- ncol(Z)
  n <- cells$n
  res <- stats::optim(c(Z, r), .circle_stress, .circle_grad, cells = cells,
                      nd = p, method = "BFGS",
                      control = list(maxit = maxit, reltol = reltol))
  Z <- matrix(res$par[seq_len(n * p)], n, p)
  r <- res$par[n * p + seq_len(n)]
  return(list(conf = Z, radius = r - min(r), stress = res$value,
              niter = unname(res$counts["gradient"]),
              converged = res$convergence == 0))
}

#' @export
print.mdspace_circle <- function(x, digits = 4, ...) {
  NextMethod()
  if (!isTRUE(x$converged)) {
    cat("Note: the optimizer did not report convergence.\n")
  }
  if (length(x$starts) > 1) {
    cat(sprintf("Starts: %d (rational + %d random); best = %s, S range %.4f to %.4f\n",
                length(x$starts), length(x$starts) - 1,
                if (x$best_start == 1) "rational" else paste("random", x$best_start - 1),
                min(x$starts), max(x$starts)))
  }
  if (isTRUE(x$threeway)) {
    cat("Two-mode three-way model (Okada & Imaizumi, 1997); S by source:\n")
    print(round(x$stress_source, digits))
    cat("Weights (w: symmetry, u: asymmetry) and aspect ratio:\n")
    wt <- x$weights
    wt$aspect <- x$aspect
    num <- vapply(wt, is.numeric, logical(1))
    wt[num] <- lapply(wt[num], round, digits)
    print(wt, row.names = FALSE)
  }
  if (length(x$path) > 1) {
    cat("S by dimensionality (rational start):\n")
    print(round(x$path, digits))
  }
  return(invisible(x))
}

#' Summarize a radius-distance (circle) model solution
#'
#' The configuration with the radii. For the two-mode three-way model the
#' weights of the sources (`w`: symmetry, `u`: asymmetry), the aspect ratio
#' of their ellipses and `S` by source are added.
#'
#' @param object An object of class `mdspace_circle`.
#' @param ... Not used.
#' @return An object of class `summary.mdspace`; for the two-mode three-way
#'   model the tables of the sources are in `extra`.
#' @references Okada, A., & Imaizumi, T. (1987). Nonmetric multidimensional
#'   scaling of asymmetric proximities. *Behaviormetrika*, 14(21), 81-96.
#'
#' Okada, A., & Imaizumi, T. (1997). Asymmetric multidimensional scaling of
#'   two-mode three-way proximities. *Journal of Classification*, 14, 195-224.
#' @export
summary.mdspace_circle <- function(object, ...) {
  out <- NextMethod()
  out$conf <- cbind(object$conf, radius = object$radius)
  if (isTRUE(object$threeway)) {
    wt <- object$weights
    wt$aspect <- unname(object$aspect)
    wt$S <- unname(object$stress_source)
    out$extra <- list("Sources (w: symmetry, u: asymmetry, aspect ratio, S)" = wt)
  }
  return(out)
}

#' Plot a radius-distance (circle) model solution
#'
#' Each object is drawn as a point with a circle of its radius. A larger
#' circle means that the object is reached more easily from the others than
#' it reaches them (for similarities).
#'
#' For the two-mode three-way model, `source` chooses what is drawn around
#' the points: the common circles of radius `r_j` (`NULL`), the ellipses of
#' one source with semiaxes `u_it r_j` (an index or a name), or the ellipses
#' of all sources overlaid in different colours (`"all"`). The axis limits
#' are the same for every choice of `source`, so that plots of different
#' sources can be compared side by side (set `par(mfrow = )` and call the
#' function once per source). In one dimension the ellipse of a source is
#' the circle of radius `u_i1 r_j`.
#'
#' @param x An object of class `mdspace_circle`.
#' @param dims Two dimensions to plot.
#' @param labels Logical; draw object labels.
#' @param source Two-mode three-way model only: `NULL` (common circles), the
#'   index or name of a source, or `"all"`.
#' @param ... Passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @references Okada, A., & Imaizumi, T. (1987). Nonmetric multidimensional
#'   scaling of asymmetric proximities. *Behaviormetrika*, 14(21), 81-96.
#'
#' Okada, A., & Imaizumi, T. (1997). Asymmetric multidimensional scaling of
#'   two-mode three-way proximities. *Journal of Classification*, 14, 195-224.
#' @export
plot.mdspace_circle <- function(x, dims = c(1, 2), labels = TRUE, source = NULL, ...) {
  conf <- x$conf
  r <- x$radius
  threeway <- isTRUE(x$threeway)
  if (!is.null(source) && !threeway) {
    stop("`source` is for the two-mode three-way model (a list of matrices).")
  }
  U <- NULL
  if (threeway) {
    U <- as.matrix(x$weights[, paste0("u", seq_len(ncol(conf))), drop = FALSE])
    rownames(U) <- x$weights$source
  }
  if (ncol(conf) == 1) {
    conf <- cbind(conf, 0)
    if (threeway) U <- cbind(U, U)
    dims <- c(1, 2)
  }
  xy <- conf[, dims, drop = FALSE]
  # どのソースを描くか(行番号にそろえる)
  idx <- integer(0)
  if (!is.null(source)) {
    if (length(source) != 1) {
      stop("`source` must be a single index, a single name, or \"all\".")
    }
    if (identical(source, "all")) {
      idx <- seq_len(nrow(U))
    } else if (is.character(source)) {
      idx <- match(source, rownames(U))
    } else {
      idx <- if (source %in% seq_len(nrow(U))) as.integer(source) else NA_integer_
    }
    if (anyNA(idx)) {
      stop("`source` does not match any source: ", source)
    }
  }
  # 軸の範囲: 二相三元は全ソースの楕円が入る大きさに固定する(ソースを替えても枠が動かない)
  if (threeway) {
    A <- outer(r, apply(U[, dims, drop = FALSE], 2, function(u) max(u, 1)))
  } else {
    A <- matrix(r, nrow(xy), 2)
  }
  lim <- range(c(xy - A, xy + A)) * 1.05
  main <- "Radius-distance model"
  if (length(idx) == 1) {
    main <- paste0(main, ": ", rownames(U)[idx])
  }
  graphics::plot(xy, asp = 1, xlim = lim, ylim = lim, pch = 16,
                 xlab = paste("Dimension", dims[1]),
                 ylab = paste("Dimension", dims[2]),
                 main = main, ...)
  graphics::abline(h = 0, v = 0, col = "grey80", lty = 3)
  th <- seq(0, 2 * pi, length.out = 181)
  # 半軸の倍率と色: 共通の円は倍率 1，ソースを選んだらその u
  if (length(idx) == 0) {
    scale <- matrix(1, 1, 2)
    cols <- "grey40"
  } else {
    scale <- U[idx, dims, drop = FALSE]
    cols <- if (length(idx) == 1) "grey40" else grDevices::hcl.colors(length(idx), "Dark 3")
  }
  for (s in seq_len(nrow(scale))) {
    for (i in seq_len(nrow(xy))) {
      if (r[i] > 0) {
        graphics::lines(xy[i, 1] + scale[s, 1] * r[i] * cos(th),
                        xy[i, 2] + scale[s, 2] * r[i] * sin(th), col = cols[s])
      }
    }
  }
  if (length(idx) > 1) {
    graphics::legend("topright", legend = rownames(U)[idx], col = cols, lty = 1,
                     bty = "n", cex = 0.8)
  }
  if (labels) {
    graphics::text(xy, labels = x$labels, pos = 3, cex = 0.8)
  }
  return(invisible(x))
}
