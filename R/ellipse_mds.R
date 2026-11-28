# 楕円モデル(ellipse-distance model): Okada (1990)，岡太・今泉の Model II と III ----------------

#' Nonmetric asymmetric MDS by the ellipse-distance (Model II) and
#' radius-distance (Model III) models
#'
#' Okada and Imaizumi (1984; see Chino & Okada, 1996, Section 3.5) wrote the
#' asymmetric proximity from object `j` to object `k` as
#' `m_jk = d_jk + alpha * c(j, k) + beta * c(k, j)`, where `d_jk` is the
#' Euclidean distance between the points `x_j` and `x_k`, and named three
#' families by how the asymmetric term `c` is restricted:
#'
#' * **Model II** (ellipse-distance model; Okada, 1990): each object `j` has an
#'   ellipse (ellipsoid, hyperellipsoid) centered at `x_j` with semiaxes
#'   `r_j1, ..., r_jp` along the dimensions, and `c(j, k)` is the radius of
#'   the ellipse of `j` in the direction of `k`, `rho_jk = d_jk v_jk` with
#'   `v_jk = 1 / sqrt(sum_t ((x_jt - x_kt) / r_jt)^2)` (Okada & Imaizumi,
#'   2024, Equations 3.10 to 3.13). When all semiaxes of `j` are equal to
#'   `r_j`, `rho_jk = r_j` and Model II reduces to Model III.
#' * **Model III** (radius-distance model; Okada & Imaizumi, 1987): `c(j, k) =
#'   r_j`, the radius of a circle (sphere, hypersphere) centered at `x_j`.
#' * **Model I**: the most general model, in which `c` depends on the pair of
#'   objects and the dimension. It is not implemented yet, because how it is
#'   restricted to be estimable is given only in Okada and Imaizumi (1984).
#'
#' The signs `alpha` and `beta` give the submodels. The defaults
#' `alpha = -1, beta = 1` are the models used by Okada and Imaizumi (1987,
#' 2024; Model IIIB): `m_jk` is the distance from the near edge of the
#' ellipse (circle) of `j` to the far edge of that of `k`, and it is assumed
#' to be monotonically decreasing in similarities (increasing in
#' dissimilarities). With `alpha = beta` the model is symmetric (Saito's
#' cluster model for Model III).
#'
#' The fit is `S = sqrt(sum((m_jk - mhat_jk)^2) / sum((m_jk - mbar)^2))`
#' with disparities `mhat_jk` from monotone regression (primary approach to
#' ties), as in [circle_mds()]. `S` is minimized over the coordinates and the
#' semiaxes (radii) by [stats::optim()] (BFGS) with the analytic gradient,
#' with the semiaxes parameterized by their logarithms so that they stay
#' positive.
#'
#' **Initial values** (Okada & Imaizumi, 2024, Section 5.4). The rational
#' start is the Torgerson configuration of the symmetrized dissimilarities
#' (with the additive constant) in `ndim_max` dimensions and the radii from
#' the skew-symmetric part, which serve as initial semiaxes for all
#' dimensions (a small constant is added so that no semiaxis is zero). The
#' solution is reduced one dimension at a time down to `ndim`: for Model II by
#' dropping the dimension with the smallest
#' `sum_j x_jt^2 + sum_j r_jt^2` (Equation 5.23), because the orientation of
#' the dimensions is determined by the model; for Model III by projecting on
#' the principal axes. `nrandom` random starts are also tried in `ndim`
#' dimensions, because the models have local minima; a typical one for Model
#' II puts most points near the origin with almost concentric ellipses
#' (Section 5.3).
#'
#' **Normalization.** The configuration is centered with the sum of squared
#' coordinates equal to `n`, and the semiaxes (radii) are rescaled with it.
#' For Model III the smallest radius is set to 0 (the radii are determined only
#' up to an additive constant); the semiaxes of Model II have no such
#' indeterminacy and are not shifted. The dimensions of Model II are not
#' rotated, because rotating them changes `S`; they are ordered by decreasing
#' variance of the coordinates. Model III is rotated to principal axes.
#'
#' **Orientation of dimensions (Model II).** If the semiaxes of two dimensions
#' are nearly equal for all objects (the points of the semiaxis configuration
#' lie near the diagonal), the ellipses are nearly circles and the orientation
#' of the two dimensions is not determined. As in Okada and Imaizumi (2024,
#' Sections 3.1.2 and 4.1.2), for each pair of dimensions the semiaxis
#' configuration is rotated by `angles` degrees counterclockwise and
#' clockwise, the coordinates are refitted with the rotated semiaxes fixed
#' (starting from the solution), and the resulting `S` is reported. If it is
#' clearly larger than the `S` of the solution, the orientation is uniquely
#' determined up to permutations and reflections of the dimensions.
#'
#' @param x A square matrix of asymmetric proximities. The diagonal is
#'   ignored; other missing cells are left out of the fit.
#' @param ndim Number of dimensions.
#' @param model `"II"` (ellipse-distance model, default) or `"III"`
#'   (radius-distance model). `"I"` is reserved and gives an error.
#' @param alpha,beta Signs of the asymmetric terms of the object itself and of
#'   the other object, each `-1`, `0` or `1` (not both 0).
#' @param type `"similarity"` (default) or `"dissimilarity"`. If `x` has a
#'   `type` attribute set by [as_proximity()], that is used.
#' @param ndim_max Dimensionality at which the rational initial configuration
#'   is built; the solutions are reduced from there down to `ndim`.
#' @param nrandom Number of random starts tried in addition to the rational
#'   start (`0` uses the rational start only).
#' @param maxit,reltol Passed to [stats::optim()].
#' @param kappa,axis_max Model II only: optional upper bounds of the aspect
#'   ratio (`max_t r_jt / min_t r_jt`) and of the size (geometric mean of the
#'   semiaxes, in units of the root mean square distance of the points from
#'   the centroid) of each ellipse. The default `Inf` leaves them unbounded as
#'   in the original method. `S` hardly determines the shape of the ellipses:
#'   without bounds some ellipses can become very long and thin or very large
#'   while `S` decreases only slightly, and with bounds the solution usually
#'   sits on the bound. A warning is given when an unbounded solution has an
#'   aspect ratio above 10 or a size above 3, and when a bound is reached.
#' @param seed Random seed for the random starts. The state of the random
#'   number generator is restored afterwards. `NULL` uses (and advances) the
#'   current state of the generator.
#' @param check_orientation Logical; for Model II with `ndim >= 2`, compute
#'   `S` after rotating the semiaxis configuration (see Details).
#' @param angles Rotation angles in degrees for `check_orientation`.
#' @return An object of class `c("mdspace_ellipse", "mdspace")` with `conf`
#'   (n x ndim), `semiaxis` (n x ndim; Model II) or `radius` (Model III),
#'   `stress` (`S`), `niter`, `mjk` and `disparity` (n x n matrices, `NA` on
#'   the diagonal and missing cells), `path` (`S` at each dimensionality from
#'   `ndim_max` to `ndim` for the rational start), `starts` (`S` of each
#'   start, the rational one first), `best_start`, `converged`, `model`,
#'   `alpha`, `beta`, `kappa`, `axis_max`, `at_bound` (whether some ellipse
#'   reached `kappa` or `axis_max`), and `orientation` (a data frame with the dimensions,
#'   angle, and `S` after counterclockwise and clockwise rotation; `NULL` when
#'   not computed).
#' @references Okada, A. (1990). A generalization of asymmetric
#'   multidimensional scaling. In M. Schader & W. Gaul (Eds.), *Knowledge,
#'   data and computer-assisted decisions* (pp. 127-138). Springer.
#'
#' Okada, A., & Imaizumi, T. (1987). Nonmetric multidimensional scaling of
#'   asymmetric proximities. *Behaviormetrika*, 14(21), 81-96.
#'
#' Okada, A., & Imaizumi, T. (2024). *Applied multidimensional scaling of
#'   asymmetric relationships*. Springer.
#'
#' 千野直仁・岡太彬訓 (1996). 非対称多次元尺度構成法とその周辺. 行動計量学,
#'   23(2), 130-152.
#' @seealso [circle_mds()] for Model IIIB with the two-mode three-way
#'   extension.
#' @examples
#' x <- matrix(c(NA, 5, 2, 1, 3,
#'               3, NA, 4, 2, 2,
#'               1, 2, NA, 6, 3,
#'               2, 1, 4, NA, 5,
#'               4, 2, 2, 3, NA), 5, byrow = TRUE)
#' fit <- ellipse_mds(x, ndim = 2, nrandom = 2)
#' fit
#' @export
ellipse_mds <- function(x, ndim = 2, model = c("II", "III", "I"), alpha = -1, beta = 1,
                        type = c("similarity", "dissimilarity"), ndim_max = ndim,
                        nrandom = 10, maxit = 2000, reltol = 1e-10, seed = 1,
                        kappa = Inf, axis_max = Inf, check_orientation = TRUE, angles = c(10, 15, 20)) {
  cl <- match.call()
  model <- match.arg(model)
  if (model == "I") {
    stop("Model I is not implemented yet: its restriction is given in Okada & Imaizumi (1984).")
  }
  if (!(alpha %in% c(-1, 0, 1)) || !(beta %in% c(-1, 0, 1)) || (alpha == 0 && beta == 0)) {
    stop("`alpha` and `beta` must each be -1, 0 or 1, and not both 0.")
  }
  if (!is.null(attr(x, "type")) && missing(type)) {
    type <- attr(x, "type")
  }
  type <- match.arg(type)
  if (!is.numeric(kappa) || length(kappa) != 1 || kappa < 1) {
    stop("`kappa` must be a number of at least 1 (`Inf` for no bound).")
  }
  if (!is.numeric(axis_max) || length(axis_max) != 1 || axis_max <= 0) {
    stop("`axis_max` must be a positive number (`Inf` for no bound).")
  }
  X <- as_proximity(x, type = type)
  n <- nrow(X)
  .stop_if_absent(X)
  if (ndim < 1 || ndim >= n) {
    stop("`ndim` must be between 1 and n - 1.")
  }
  if (ndim_max < ndim || ndim_max >= n) {
    stop("`ndim_max` must be between `ndim` and n - 1.")
  }

  # 当てはめに使うセル(対角と欠測を除く) ---------------------------------------
  use <- row(X) != col(X) & !is.na(X)
  cells <- list(s = X[use], j = row(X)[use], k = col(X)[use], n = n,
                decreasing = type == "similarity", model = model,
                alpha = alpha, beta = beta, kappa = kappa, axis_max = axis_max)
  if (length(cells$s) < 3) {
    stop("Too few observed off-diagonal cells.")
  }

  # 合理的初期布置: 最大次元で作り，1次元ずつ下ろす -------------------------------------
  init <- .circle_init(X, ndim_max, type)
  Z <- init$conf
  if (model == "II") {
    # 円モデルの初期半径を全次元の半軸に(5.4.4)。半軸が 0 だと動かないので小さな値を足す
    s0 <- sqrt(mean(rowSums(sweep(Z, 2, colMeans(Z))^2)))
    A <- matrix(init$radius + 0.1 * s0, n, ndim_max)
  } else {
    A <- init$radius
  }
  path <- numeric(0)
  fit <- NULL
  for (p in seq(ndim_max, ndim)) {
    if (ncol(Z) > p) {
      if (model == "II") {
        # 式(5.23): 座標と半軸の二乗和の大きい p 次元を残す(回転はしない)
        Zc <- sweep(Z, 2, colMeans(Z))
        keep <- sort(order(colSums(Zc^2) + colSums(A^2), decreasing = TRUE)[seq_len(p)])
        Z <- Zc[, keep, drop = FALSE]
        A <- A[, keep, drop = FALSE]
      } else {
        Z <- .principal_axes(Z, p)
      }
    }
    fit <- .ell_optim(Z, A, cells, maxit, reltol)
    Z <- fit$conf
    A <- fit$axis
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
      A0 <- if (model == "II") matrix(stats::runif(n * ndim, 0.1, 0.6), n, ndim) else stats::runif(n, 0, 0.5)
      f <- .ell_optim(Z0, A0, cells, maxit, reltol)
      starts <- c(starts, f$stress)
      if (f$stress < best$stress) {
        best <- f
        best_start <- i + 1L
      }
    }
  }

  # 正規化して返す ---------------------------------------------------------------
  nr <- .ell_normalize(best$conf, best$axis, model)
  conf <- nr$conf
  rownames(conf) <- rownames(X)
  colnames(conf) <- paste0("D", seq_len(ndim))
  axis <- nr$axis
  if (model == "II") {
    dimnames(axis) <- list(rownames(X), paste0("D", seq_len(ndim)))
  } else {
    names(axis) <- rownames(X)
  }
  M <- D <- matrix(NA_real_, n, n, dimnames = dimnames(X))
  M[use] <- .ell_parts(.ell_par(conf, axis, model), cells, ndim)$m
  D[use] <- .monotone(M[use], cells$s, cells$decreasing)

  # 次元の向きの確かめ(Model II・2次元以上) --------------------------------------------
  # 原典どおり既定では縛らない。そのかわり極端な楕円が出たら警告する
  # (目安: 縦横比 10 倍・大きさ 3。本の表4.4 の解は縦横比 3.3・大きさ 1.18)
  if (model == "II") {
    ratio <- apply(axis, 1, max) / apply(axis, 1, min)
    gsize <- exp(rowMeans(log(axis)))
    if (!is.finite(kappa) && max(ratio) > 10) {
      warning(sprintf("Some ellipses are very long and thin (largest aspect ratio %.3g, e.g. %s); the shape of the ellipses is weakly determined by S. Consider `kappa`.",
                      max(ratio), rownames(X)[which.max(ratio)]))
    }
    if (!is.finite(axis_max) && max(gsize) > 3) {
      warning(sprintf("Some ellipses are very large relative to the configuration (largest size %.3g, e.g. %s); the shape of the ellipses is weakly determined by S. Consider `axis_max`.",
                      max(gsize), rownames(X)[which.max(gsize)]))
    }
  }
  if (isTRUE(best$at_kappa)) {
    warning(sprintf("The aspect ratio of some ellipses reached the bound `kappa` = %g; the data may favour collapsed ellipses.", kappa))
  }
  if (isTRUE(best$at_size)) {
    warning(sprintf("The size of some ellipses reached the bound `axis_max` = %g; the data may favour very large ellipses.", axis_max))
  }
  orientation <- NULL
  if (model == "II" && ndim >= 2 && isTRUE(check_orientation)) {
    orientation <- .ell_orientation(conf, axis, cells, angles, maxit, reltol)
  }

  out <- .new_mdspace("ellipse", conf = conf, stress = best$stress,
                      niter = best$niter, data = X, call = cl, ndim = ndim,
                      mjk = M, disparity = D, type = type, path = path,
                      starts = starts, best_start = best_start,
                      converged = best$converged, variant = model,
                      alpha = alpha, beta = beta, orientation = orientation,
                      kappa = if (model == "II") kappa else NULL,
                      axis_max = if (model == "II") axis_max else NULL,
                      at_bound = c(kappa = isTRUE(best$at_kappa), axis_max = isTRUE(best$at_size)))
  if (model == "II") {
    out$semiaxis <- axis
  } else {
    out$radius <- axis
  }
  return(out)
}

# 母数ベクトルの組み立て(内部) ----------------------------------------------------------
# Model II は半軸を対数で持つ(正に保つため)，Model III は半径をそのまま持つ
.ell_par <- function(Z, A, model) {
  if (model == "II") {
    return(c(Z, log(A)))
  }
  return(c(Z, A))
}

# m_jk とその部品(内部) ----------------------------------------------------------------
# rho_jk = 楕円 j の k の方向の半径 = d_jk v_jk，v_jk = 1 / sqrt(sum_t (dif_t / r_jt)^2)(dif = x_j - x_k)
.ell_parts <- function(par, cells, nd) {
  p <- nd
  n <- cells$n
  Z <- matrix(par[seq_len(n * p)], n, p)
  dif <- Z[cells$j, , drop = FALSE] - Z[cells$k, , drop = FALSE]
  d <- sqrt(rowSums(dif^2))
  if (cells$model == "II") {
    # 対数の半軸が極端に振れて 0 や Inf にならないよう，exp(±30) で止める
    R <- matrix(exp(pmin(pmax(par[n * p + seq_len(n * p)], -30), 30)), n, p)
    Rj <- R[cells$j, , drop = FALSE]
    Rk <- R[cells$k, , drop = FALSE]
    # 重なった点(d = 0)では向きが決まらないので，半径は 0 とする
    ok <- d > 1e-12
    vj <- vk <- numeric(length(d))
    vj[ok] <- 1 / sqrt(rowSums((dif[ok, , drop = FALSE] / Rj[ok, , drop = FALSE])^2))
    vk[ok] <- 1 / sqrt(rowSums((dif[ok, , drop = FALSE] / Rk[ok, , drop = FALSE])^2))
    # 式(3.10): 方向の半径は d_jk v_jk
    rj <- d * vj
    rk <- d * vk
  } else {
    r <- par[n * p + seq_len(n)]
    Rj <- Rk <- vj <- vk <- NULL
    rj <- r[cells$j]
    rk <- r[cells$k]
  }
  m <- d + cells$alpha * rj + cells$beta * rk
  return(list(dif = dif, d = d, vj = vj, vk = vk, rj = rj, rk = rk, Rj = Rj, Rk = Rk, m = m))
}

.ell_stress <- function(par, cells, nd) {
  q <- .ell_parts(par, cells, nd)
  mh <- .monotone(q$m, cells$s, cells$decreasing)
  return(sqrt(sum((q$m - mh)^2) / sum((q$m - mean(q$m))^2)))
}

# S の勾配(ディスパリティを固定した解析的な勾配)(内部) -----------------------------------------
.ell_grad <- function(par, cells, nd) {
  n <- cells$n
  p <- nd
  q <- .ell_parts(par, cells, p)
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
  al <- cells$alpha
  be <- cells$beta
  idx <- factor(c(cells$j, cells$k), levels = seq_len(n))
  wd <- ifelse(q$d > 1e-12, 1 / q$d, 0)
  GZ <- matrix(0, n, p)
  if (cells$model == "II") {
    GA <- matrix(0, n, p)
    # rho = d v，v = Q^(-1/2)，Q = sum_t (dif_t / r_t)^2 として
    # d rho / d dif_t = dif_t v / d - d v^3 dif_t / r_t^2(rho_kj も dif の向きによらず同じ形)
    # d rho / d log r_t = d v^3 dif_t^2 / r_t^2
    dvj3 <- q$d * q$vj^3
    dvk3 <- q$d * q$vk^3
    for (t in seq_len(p)) {
      dt <- q$dif[, t]
      # dm/d(dif_t): 距離の項と，2つの方向の半径の項
      dm <- dt * wd +
        al * (dt * q$vj * wd - dvj3 * dt / q$Rj[, t]^2) +
        be * (dt * q$vk * wd - dvk3 * dt / q$Rk[, t]^2)
      v <- g * dm
      GZ[, t] <- as.vector(tapply(c(v, -v), idx, sum, default = 0))
      # dm/d(log r_jt) と dm/d(log r_kt)
      aj <- g * al * dvj3 * dt^2 / q$Rj[, t]^2
      ak <- g * be * dvk3 * dt^2 / q$Rk[, t]^2
      GA[, t] <- as.vector(tapply(c(aj, ak), idx, sum, default = 0))
    }
    return(c(GZ, GA))
  }
  for (t in seq_len(p)) {
    v <- g * q$dif[, t] * wd
    GZ[, t] <- as.vector(tapply(c(v, -v), idx, sum, default = 0))
  }
  gr <- as.vector(tapply(c(al * g, be * g), idx, sum, default = 0))
  return(c(GZ, gr))
}

# 1つの初期値からの最小化(内部) -----------------------------------------------------
# Model II は log r_jt = c_j + e_jt と分け，形 e_jt だけを [-L, L](L = log(kappa) / 2)に収めて
# L-BFGS-B で解く＝各対象の縦横比 max_t r_jt / min_t r_jt <= kappa(楕円が線分につぶれる退化を止める)。
# e_jt は平均 0 に直して使うので，c_j は幾何平均の対数(axis_max の上限はここに置く)
.ell_optim <- function(Z, A, cells, maxit, reltol) {
  p <- ncol(Z)
  n <- cells$n
  if (cells$model == "III") {
    res <- stats::optim(.ell_par(Z, A, "III"), .ell_stress, .ell_grad,
                        cells = cells, nd = p, method = "BFGS",
                        control = list(maxit = maxit, reltol = reltol))
    return(list(conf = matrix(res$par[seq_len(n * p)], n, p),
                axis = res$par[n * p + seq_len(n)], stress = res$value,
                niter = unname(res$counts["gradient"]),
                converged = res$convergence == 0, at_kappa = FALSE, at_size = FALSE))
  }
  L <- log(cells$kappa) / 2
  U <- log(cells$axis_max)
  nz <- n * p
  # 座標の大きさは目的関数の中で固定する(重心からの二乗平均平方根 = 1，二乗和 = n)。
  # S は座標と半軸を同じ比で伸ばしても変わらないので，こうしないと半軸の大きさの上限が効かない
  scl <- function(z) {
    Zc <- sweep(matrix(z, n, p), 2, colMeans(matrix(z, n, p)))
    sqrt(n / sum(Zc^2))
  }
  # 初期の半軸も，正規化した座標の単位に直す
  eta0 <- log(A) + log(scl(Z))
  c0 <- pmin(rowMeans(eta0), U)
  E0 <- pmin(pmax(eta0 - rowMeans(eta0), -L), L)
  to_full <- function(th) {
    z <- th[seq_len(nz)]
    # 形 e_jt は対象ごとに平均 0 に直して足す＝c_j がちょうど半軸の幾何平均の対数になる
    E <- matrix(th[nz + n + seq_len(nz)], n, p)
    c(scl(z) * z, th[nz + seq_len(n)] + E - rowMeans(E))
  }
  fn <- function(th) .ell_stress(to_full(th), cells, p)
  gr <- function(th) {
    z <- th[seq_len(nz)]
    a <- scl(z)
    g <- .ell_grad(to_full(th), cells, p)
    GZ <- g[seq_len(nz)]
    Zc <- c(sweep(matrix(z, n, p), 2, colMeans(matrix(z, n, p))))
    # Zs = a(Z) Z の連鎖律: a (G - (sum(G Zc) / sum(Zc^2)) Zc)
    GZ <- a * (GZ - sum(GZ * Zc) / sum(Zc^2) * Zc)
    GA <- matrix(g[nz + seq_len(nz)], n, p)
    c(GZ, rowSums(GA), GA - rowMeans(GA))
  }
  th0 <- c(Z, c0, E0)
  if (is.finite(L) || is.finite(U)) {
    lower <- c(rep(-Inf, nz + n), rep(-L, nz))
    upper <- c(rep(Inf, nz), rep(U, n), rep(L, nz))
    res <- stats::optim(th0, fn, gr, method = "L-BFGS-B", lower = lower, upper = upper,
                        control = list(maxit = maxit, factr = reltol / .Machine$double.eps))
  } else {
    res <- stats::optim(th0, fn, gr, method = "BFGS",
                        control = list(maxit = maxit, reltol = reltol))
  }
  full <- to_full(res$par)
  A <- matrix(exp(pmin(pmax(full[nz + seq_len(nz)], -30), 30)), n, p)
  E <- matrix(res$par[nz + n + seq_len(nz)], n, p)
  cc <- res$par[nz + seq_len(n)]
  # 上限に届いたか: 縦横比(ある対象の形の幅が 2L)と大きさ(c_j が U)
  at_kappa <- is.finite(L) && p > 1 && any(apply(E, 1, function(e) max(e) - min(e)) > 2 * L - 1e-6)
  at_size <- is.finite(U) && any(cc > U - 1e-6)
  return(list(conf = matrix(full[seq_len(nz)], n, p), axis = A, stress = res$value,
              niter = unname(res$counts["gradient"]),
              converged = res$convergence == 0, at_kappa = at_kappa, at_size = at_size))
}

# 正規化(内部) -------------------------------------------------------------------------
# 重心を原点・座標の二乗和を n・半軸(半径)も同じ比で伸縮。
# Model II は回転せず，座標の分散の大きい順に次元を並べる。Model III は主軸に回し，最小の半径を 0 に
.ell_normalize <- function(Z, A, model) {
  Z <- sweep(Z, 2, colMeans(Z))
  a <- sqrt(nrow(Z) / sum(Z^2))
  Z <- Z * a
  A <- A * a
  if (model == "II") {
    o <- order(colSums(Z^2), decreasing = TRUE)
    return(list(conf = Z[, o, drop = FALSE], axis = A[, o, drop = FALSE]))
  }
  Z <- Z %*% svd(Z)$v
  return(list(conf = Z, axis = A - min(A)))
}

# 次元の向きの確かめ(内部) ------------------------------------------------------------
# 半軸の布置の2次元 (t, s) を反時計回り・時計回りに回し，回した半軸を固定して座標だけを最小化する。
# 半軸は2乗でしか効かないので，回して負になった値は絶対値にする
.ell_orientation <- function(conf, axis, cells, angles, maxit, reltol) {
  p <- ncol(conf)
  n <- nrow(conf)
  # 半軸を固定した S と勾配(座標の部分だけ)
  fz <- function(z, la) .ell_stress(c(z, la), cells, p)
  gz <- function(z, la) .ell_grad(c(z, la), cells, p)[seq_len(n * p)]
  pairs <- utils::combn(p, 2)
  rows <- list()
  for (k in seq_len(ncol(pairs))) {
    ts <- pairs[, k]
    for (ang in angles) {
      sv <- vapply(c(1, -1), function(sgn) {
        th <- sgn * ang * pi / 180
        R <- axis
        Rot <- matrix(c(cos(th), sin(th), -sin(th), cos(th)), 2)
        # 行ベクトル (r_t, r_s) を反時計回りに th 回す
        R[, ts] <- abs(axis[, ts] %*% t(Rot))
        res <- stats::optim(c(conf), fz, gz, la = c(log(R)), method = "BFGS",
                            control = list(maxit = maxit, reltol = reltol))
        res$value
      }, numeric(1))
      rows[[length(rows) + 1]] <- data.frame(dims = paste(ts, collapse = "-"), angle = ang,
                                             counterclockwise = sv[1], clockwise = sv[2])
    }
  }
  return(do.call(rbind, rows))
}

#' @export
print.mdspace_ellipse <- function(x, digits = 4, ...) {
  NextMethod()
  lab <- if (x$variant == "II") "ellipse-distance" else "radius-distance"
  cat(sprintf("Model %s (%s), alpha = %d, beta = %d\n", x$variant, lab, x$alpha, x$beta))
  if (!isTRUE(x$converged)) {
    cat("Note: the optimizer did not report convergence.\n")
  }
  if (length(x$starts) > 1) {
    cat(sprintf("Starts: %d (rational + %d random); best = %s, S range %.4f to %.4f\n",
                length(x$starts), length(x$starts) - 1,
                if (x$best_start == 1) "rational" else paste("random", x$best_start - 1),
                min(x$starts), max(x$starts)))
  }
  if (length(x$path) > 1) {
    cat("S by dimensionality (rational start):\n")
    print(round(x$path, digits))
  }
  if (!is.null(x$orientation)) {
    cat(sprintf("S after rotating the semiaxis configuration (solution S = %.4f):\n", x$stress))
    o <- x$orientation
    o$counterclockwise <- round(o$counterclockwise, digits)
    o$clockwise <- round(o$clockwise, digits)
    print(o, row.names = FALSE)
  }
  return(invisible(x))
}

#' @export
summary.mdspace_ellipse <- function(object, ...) {
  out <- NextMethod()
  if (object$variant == "II") {
    sa <- object$semiaxis
    colnames(sa) <- paste0("r", seq_len(ncol(sa)))
    out$conf <- cbind(object$conf, sa)
  } else {
    out$conf <- cbind(object$conf, radius = object$radius)
  }
  return(out)
}

#' Plot an ellipse-distance or radius-distance model solution
#'
#' `which = "conf"` draws each object as a point with its ellipse (Model II)
#' or circle (Model III) in two dimensions. `which = "semiaxis"` draws the
#' semiaxis configuration of Model II with the diagonal; when the points lie
#' near the diagonal the ellipses are nearly circles and the orientation of
#' the two dimensions is not determined (Okada & Imaizumi, 2024, Fig. 3.4).
#'
#' @param x An object of class `mdspace_ellipse`.
#' @param which `"conf"` or `"semiaxis"`.
#' @param dims Two dimensions to plot.
#' @param labels Logical; draw object labels.
#' @param ... Passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
plot.mdspace_ellipse <- function(x, which = c("conf", "semiaxis"), dims = c(1, 2),
                                 labels = TRUE, ...) {
  which <- match.arg(which)
  conf <- x$conf
  n <- nrow(conf)
  if (x$variant == "II") {
    sa <- x$semiaxis
  } else {
    sa <- matrix(x$radius, n, ncol(conf))
  }
  if (ncol(conf) == 1) {
    conf <- cbind(conf, 0)
    sa <- cbind(sa, sa)
    dims <- c(1, 2)
  }
  if (which == "semiaxis") {
    if (x$variant != "II") {
      stop("The semiaxis configuration is for Model II.")
    }
    xy <- sa[, dims, drop = FALSE]
    lim <- c(0, max(xy) * 1.1)
    graphics::plot(xy, xlim = lim, ylim = lim, asp = 1, pch = 16,
                   xlab = paste("Semiaxis along dimension", dims[1]),
                   ylab = paste("Semiaxis along dimension", dims[2]),
                   main = "Semiaxis configuration", ...)
    graphics::abline(0, 1, col = "grey60", lty = 2)
    if (labels) {
      graphics::text(xy, labels = x$labels, pos = 3, cex = 0.8)
    }
    return(invisible(x))
  }
  xy <- conf[, dims, drop = FALSE]
  ab <- sa[, dims, drop = FALSE]
  lim <- range(c(xy - ab, xy + ab)) * 1.05
  main <- if (x$variant == "II") "Ellipse-distance model" else "Radius-distance model"
  graphics::plot(xy, asp = 1, xlim = lim, ylim = lim, pch = 16,
                 xlab = paste("Dimension", dims[1]),
                 ylab = paste("Dimension", dims[2]), main = main, ...)
  graphics::abline(h = 0, v = 0, col = "grey80", lty = 3)
  th <- seq(0, 2 * pi, length.out = 181)
  for (i in seq_len(n)) {
    if (all(ab[i, ] > 0)) {
      graphics::lines(xy[i, 1] + ab[i, 1] * cos(th), xy[i, 2] + ab[i, 2] * sin(th),
                      col = "grey40")
    }
  }
  if (labels) {
    graphics::text(xy, labels = x$labels, pos = 3, cex = 0.8)
  }
  return(invisible(x))
}
