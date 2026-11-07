# GSTATIS(二相三元の HFM; Chino, Grorud & Yoshino, 1996) --------------------------

#' GSTATIS: STATIS for a set of asymmetric matrices via Hermitian forms
#'
#' A set of `q` square asymmetric matrices on the same `n` objects (e.g. the
#' same group rated at `q` occasions) is analyzed in three steps, following
#' the revised version of GSTATIS (Chino, Grorud, & Yoshino).
#'
#' 1. **Inter-structure.** Each matrix is turned into its Hermitian form
#'    `H_i = S_i + iK_i`. The generalized RV coefficient between sources is
#'    `GRV(H_i, H_i') = tr(H_i H_i') / sqrt(tr(H_i^2) tr(H_i'^2))`, and the
#'    `q x q` GRV matrix is eigen-decomposed.
#' 2. **Intra-structure.** The `a`-th eigenvector `beta_a` of the GRV matrix
#'    gives the weights of the compromise `Xi_a = sum_i beta_ia H_i`, which is
#'    analyzed by HFM (see [hfm_mds()]). The sign of an eigenvector is
#'    arbitrary: `beta_1` is signed so that its sum is positive, and each
#'    later `beta_a` so that its largest element in absolute value is
#'    positive. Flipping `beta_a` negates `Xi_a`, which swaps positive and
#'    negative planes and so turns the trajectories by 90 degrees.
#' 3. **Trajectories (revised step 3).** The `j`-th row (or column) vector
#'    `h_j` of each `H_i` is projected onto the compromise configuration
#'    `O_p = U |Psi|^{1/2}`: the image is `Psi^{-1/2} O_p* h_j`, where a
#'    negative eigenvalue `psi` is read as `sqrt(psi) = i sqrt(|psi|)`. With
#'    unit eigenvectors `U` this is `D U* h_j` with `D = 1` for a positive and
#'    `-i` for a negative eigenvalue. The distance between two sources for
#'    object `j`, `||Psi^{-1/2} O_p* (h_j^(i) - h_j^(i'))||`, is returned in
#'    `lpnorm`.
#'
#' The results reproduce the FORTRAN program by N. Chino (output of
#' 1999-12-30): the GRV matrix, its eigenvalues and eigenvectors, the
#' eigenvalues of the compromises, and all the trajectories (the program
#' prints each object's trajectory divided by its largest norm over sources;
#' see `plot(..., scale = "max")`).
#'
#' @param x A list of square matrices with the same labels, or an
#'   `n x n x T` array.
#' @param ndim Number of complex dimensions to keep.
#' @param ncomp Number of compromises (points of view) to build. Use
#'   `ncomp = NULL` with `threshold` to keep every RV component whose
#'   eigenvalue exceeds `threshold`, as in the 2001 program.
#' @param threshold Eigenvalue threshold used when `ncomp = NULL`.
#' @param scale `"none"` (default, as in Chino's program) or `"norm"` (each
#'   `H_i` divided by its Frobenius norm before step 2, the usual STATIS
#'   normalization).
#' @param project `"row"` (default, as in Chino's program) or `"column"`:
#'   which vectors of `H_i` are projected in step 3.
#' @param missing,impute_ndim,tol,maxit Handling of missing values in each
#'   source, as in [hfm_mds()]. An object absent at a source (whole row and
#'   column missing) is an error.
#' @param rotate Rotation of each complex plane; see [hfm_mds()].
#' @return An object of class `c("mdspace_gstatis", "mdspace")` with
#'   `conf`, `coord`, `eigen`, `gof`, `unexplained` (as in [hfm_mds()], for the
#'   first compromise), `rv` (RV matrix), `rv_eigen` (its eigenvalues and
#'   vectors), `weights` (`q x ncomp`), `compromise` (list, one per point of
#'   view, each with the HFM of `Xi_a`, its `trajectory`, `lpnorm`, and
#'   `distance`, the Frobenius distance `||H_i - Xi_a||` for each source),
#'   `trajectory` and `lpnorm` (those of the first point of view; complex
#'   array `n x ndim x q` and array `q x q x n`), and `occasion` (HFM of each
#'   source, see [hfm_mds()]).
#' @references
#' Chino, N., Grorud, A., & Yoshino, R. (1996). A complex analysis for
#'   two-mode three-way asymmetric relational data. *Abstracts of the talks at
#'   the 5th Conference of the International Federation of Classification
#'   Societies*, Vol. 2, 83-86.
#'
#' Chino, N., Grorud, A., & Yoshino, R. A revised version of GSTATIS.
#'   *Proceedings of the Annual Meeting of the Behaviormetric Society of
#'   Japan*, 409-412.
#'
#' Lavit, C., Escoufier, Y., Sabatier, R., & Traissac, P. (1994). The ACT
#'   (STATIS method). *Computational Statistics & Data Analysis*, 18, 97-119.
#' @examples
#' set.seed(1)
#' base <- matrix(sample(1:7, 25, TRUE), 5, 5, dimnames = list(LETTERS[1:5], LETTERS[1:5]))
#' xs <- lapply(1:3, function(t) base + matrix(sample(-1:1, 25, TRUE), 5, 5))
#' fit <- gstatis_mds(xs, ndim = 2)
#' fit
#' plot(fit, type = "trajectory", dim = 1)
#' @export
gstatis_mds <- function(x, ndim = 2, ncomp = 1, threshold = 1,
                        scale = c("none", "norm"),
                        rotate = c("std", "varimax", "none"),
                        project = c("row", "column"),
                        missing = c("none", "impute"), impute_ndim = ndim,
                        tol = 1e-4, maxit = 1000) {
  cl <- match.call()
  missing <- match.arg(missing)
  scale <- match.arg(scale)
  project <- match.arg(project)
  rotate <- match.arg(rotate)

  # 入力を正方行列のリストにそろえる -------------------------------------------------
  if (is.array(x) && length(dim(x)) == 3) {
    x <- lapply(seq_len(dim(x)[3]), function(t) x[, , t])
  }
  if (!is.list(x) || length(x) < 2) {
    stop("`x` must be a list of at least two square matrices or an n x n x T array.")
  }
  X <- lapply(x, as_proximity, type = "similarity")
  n <- nrow(X[[1]])
  if (any(vapply(X, nrow, integer(1)) != n)) {
    stop("All matrices must have the same size.")
  }
  labs <- rownames(X[[1]])
  X <- lapply(X, function(m) {
    dimnames(m) <- list(labs, labs)
    return(m)
  })
  occ <- names(x)
  if (is.null(occ)) occ <- paste0("T", seq_along(X))
  nT <- length(X)

  # 欠測の扱い: 時点ごとに MINERS 法で補完する ----------------------------------------
  imputed <- NULL
  if (any(vapply(X, anyNA, logical(1)))) {
    if (missing == "none") {
      stop("`x` has missing values; use missing = \"impute\" (see missing_pattern()).")
    }
    imputed <- list()
    for (t in seq_len(nT)) {
      if (anyNA(X[[t]])) {
        imp <- .impute_miners(X[[t]], p = impute_ndim, tol = tol, maxit = maxit, label = occ[t])
        X[[t]] <- imp$data
        imputed[[occ[t]]] <- cbind(source = occ[t], imp$imputed)
      }
    }
    imputed <- do.call(rbind, imputed)
    rownames(imputed) <- NULL
  }
  if (ndim < 1 || ndim > n) {
    stop("`ndim` must be between 1 and n.")
  }

  # 第1段: 時点間の一般化 RV 係数 ----------------------------------------------------
  H <- lapply(X, .hermitian)
  if (scale == "norm") {
    H <- lapply(H, function(h) h / sqrt(sum(Mod(h)^2)))
  }
  ip <- matrix(0, nT, nT)
  for (s in seq_len(nT)) {
    for (t in seq_len(nT)) {
      # tr(H_s H_t) = sum(H_s * Conj(H_t)) (エルミート行列なので実数)
      ip[s, t] <- Re(sum(H[[s]] * Conj(H[[t]])))
    }
  }
  rv <- ip / sqrt(outer(diag(ip), diag(ip)))
  dimnames(rv) <- list(occ, occ)
  rv_eig <- eigen(rv, symmetric = TRUE)
  V <- rv_eig$vectors
  # 固有ベクトルの符号は不定なので決めておく: 第1は和が正，第2以降は絶対値が最大の要素が正
  if (sum(V[, 1]) < 0) {
    V[, 1] <- -V[, 1]
  }
  for (k in seq_len(nT)[-1]) {
    if (V[which.max(abs(V[, k])), k] < 0) {
      V[, k] <- -V[, k]
    }
  }
  # 妥協の重みは GRV の固有ベクトルそのもの(論文の式3・千野の FORTRAN 版と同じ)
  dimnames(V) <- list(occ, paste0("C", seq_len(nT)))

  if (is.null(ncomp)) {
    ncomp <- max(1, sum(rv_eig$values > threshold))
  }
  ncomp <- min(ncomp, nT)

  # 第2段: 視点ごとの妥協行列を HFM にかける ------------------------------------------
  # 第3段: 各時点の行(列)ベクトルを妥協の配置へ射影する(改訂版の式5) -------------------
  keep <- seq_len(ndim)
  comp <- vector("list", ncomp)
  for (k in seq_len(ncomp)) {
    Xi <- Reduce(`+`, Map(function(h, w) h * w, H, V[, k]))
    eg <- .hfm_eigen(Xi, rotate = rotate)
    lam <- eg$values
    coord <- eg$vectors[, keep, drop = FALSE] %*% diag(sqrt(abs(lam[keep])), ndim)
    dimnames(coord) <- list(labs, paste0("D", keep))
    gof <- lam^2 / sum(lam^2)
    # Psi^{-1/2} O* h = D U* h(負の固有値は sqrt(psi) = i sqrt|psi| と読むので D = -i)
    U <- eg$vectors[, keep, drop = FALSE]
    D <- ifelse(lam[keep] < 0, -1i, 1 + 0i)
    traj <- array(0i, dim = c(n, ndim, nT),
                  dimnames = list(labs, paste0("D", keep), occ))
    for (t in seq_len(nT)) {
      # 行ベクトル h_j の像: sum_l conj(u_l) H[j, l]。列なら sum_l conj(u_l) H[l, j]
      Ht <- if (project == "row") H[[t]] else t(H[[t]])
      traj[, , t] <- (Ht %*% Conj(U)) %*% diag(D, ndim)
    }
    # 対象ごとの時点間の距離(Lp ノルム)
    lp <- array(0, dim = c(nT, nT, n), dimnames = list(occ, occ, labs))
    for (j in seq_len(n)) {
      Z <- matrix(traj[j, , ], nrow = ndim)
      for (s1 in seq_len(nT)) {
        for (s2 in seq_len(nT)) {
          lp[s1, s2, j] <- sqrt(sum(Mod(Z[, s1] - Z[, s2])^2))
        }
      }
    }
    # 各時点のエルミート行列と妥協行列のフロベニウス距離
    dist_xi <- stats::setNames(vapply(H, function(h) sqrt(sum(Mod(h - Xi)^2)), numeric(1)), occ)
    comp[[k]] <- list(matrix = Xi, weights = V[, k], vectors = eg$vectors,
                      eigen = lam, coord = coord, gof = gof,
                      unexplained = 1 - sum(gof[keep]),
                      trajectory = traj, lpnorm = lp, distance = dist_xi)
  }
  names(comp) <- paste0("C", seq_len(ncomp))
  c1 <- comp[[1]]

  # 準備段階: 各時点の HFM ---------------------------------------------------------
  occasion <- lapply(X, hfm_mds, ndim = ndim, rotate = rotate)
  names(occasion) <- occ

  out <- .new_mdspace("gstatis", conf = .complex_to_conf(c1$coord), data = X,
                      call = cl, ndim = ndim, coord = c1$coord,
                      eigen = c1$eigen, gof = c1$gof,
                      unexplained = c1$unexplained,
                      rv = rv, rv_eigen = list(values = rv_eig$values, vectors = V),
                      weights = V[, seq_len(ncomp), drop = FALSE],
                      compromise = comp, trajectory = c1$trajectory,
                      lpnorm = c1$lpnorm, occasion = occasion,
                      scale = scale, project = project)
  if (!is.null(imputed)) {
    out$imputed <- imputed
  }
  return(out)
}

#' @export
print.mdspace_gstatis <- function(x, digits = 4, ...) {
  NextMethod()
  cat("Occasions:", dim(x$trajectory)[3],
      " Points of view:", length(x$compromise), "\n")
  cat("\nRV eigenvalues (interstructure):\n")
  ev <- x$rv_eigen$values
  print(round(rbind(eigenvalue = ev, share = ev / sum(ev)), digits))
  cat("\nWeights of the first compromise (GRV eigenvector):\n")
  print(round(x$weights[, 1], digits))
  k <- ncol(x$coord)
  cat("\nFirst compromise, Unexplained SS (1 - cumulative GOF):",
      format(round(x$unexplained, digits), nsmall = digits), "\n")
  tab <- rbind(eigenvalue = x$eigen[seq_len(k)], gof = x$gof[seq_len(k)],
               cumulative = cumsum(x$gof[seq_len(k)]))
  colnames(tab) <- colnames(x$coord)
  print(round(tab, digits))
  return(invisible(x))
}

#' Plot a GSTATIS solution
#'
#' `type = "interstructure"` plots the sources on the first two eigenvectors
#' of the GRV matrix. `type = "compromise"` plots one complex plane of a
#' compromise like [plot.mdspace_hfm()]. `type = "trajectory"` draws, for each
#' object, its images over the sources joined by lines (step 3); the filled
#' point is the image of the compromise row, i.e. the weighted mean of the
#' images with the GRV weights.
#'
#' @param x An object of class `mdspace_gstatis`.
#' @param type `"trajectory"`, `"compromise"`, or `"interstructure"`.
#' @param dim Which complex dimension to plot.
#' @param view Which compromise (point of view).
#' @param objects Objects to draw in `type = "trajectory"` (names or
#'   indices); all by default.
#' @param scale `"none"` or `"max"`: with `"max"` each object's trajectory is
#'   divided by its largest norm over sources (over all kept dimensions), as
#'   printed by Chino's program.
#' @param ... Passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
plot.mdspace_gstatis <- function(x, type = c("trajectory", "compromise", "interstructure"),
                                 dim = 1, view = 1, objects = NULL,
                                 scale = c("none", "max"), ...) {
  type <- match.arg(type)
  scale <- match.arg(scale)
  if (type == "interstructure") {
    V <- x$rv_eigen$vectors
    ev <- x$rv_eigen$values
    xy <- V[, 1:2, drop = FALSE] %*% diag(sqrt(pmax(ev[1:2], 0)), 2)
    lim <- c(-1, 1) * max(abs(xy)) * 1.2
    graphics::plot(xy, asp = 1, xlim = lim, ylim = lim, pch = 16,
                   xlab = "Component 1", ylab = "Component 2",
                   main = "GSTATIS inter-structure", ...)
    graphics::abline(h = 0, v = 0, col = "grey80", lty = 3)
    graphics::text(xy, labels = rownames(V), pos = 3, cex = 0.8)
    return(invisible(x))
  }
  cp <- x$compromise[[view]]
  if (type == "compromise") {
    z <- cp$coord[, dim]
    lim <- c(-1, 1) * max(Mod(z)) * 1.2
    graphics::plot(Re(z), Im(z), asp = 1, xlim = lim, ylim = lim, type = "n",
                   xlab = "Real", ylab = "Imaginary",
                   main = sprintf("GSTATIS compromise %d, dimension %d (eigenvalue %.3f)",
                                  view, dim, cp$eigen[dim]), ...)
    graphics::abline(h = 0, v = 0, col = "grey80", lty = 3)
    ok <- Mod(z) > 1e-3 * max(Mod(z))
    graphics::arrows(0, 0, Re(z[ok]), Im(z[ok]), length = 0.08, col = "grey40")
    graphics::text(Re(z), Im(z), labels = x$labels, pos = 3, cex = 0.8)
    return(invisible(x))
  }
  # 軌跡: 時点ごとの像を線でつなぎ，妥協の行の像(重みつき平均)を黒丸で示す
  arr <- cp$trajectory
  if (scale == "max") {
    for (j in seq_len(dim(arr)[1])) {
      m <- max(sqrt(colSums(Mod(matrix(arr[j, , ], nrow = dim(arr)[2]))^2)))
      if (m > 0) arr[j, , ] <- arr[j, , ] / m
    }
  }
  tr <- matrix(arr[, dim, ], nrow = dim(arr)[1], dimnames = list(x$labels, dimnames(arr)[[3]]))
  w <- cp$weights
  z0 <- as.vector(tr %*% w) / sum(w)
  sel <- if (is.null(objects)) seq_len(nrow(tr)) else objects
  tr <- tr[sel, , drop = FALSE]
  z0 <- stats::setNames(z0, x$labels)[sel]
  allz <- c(tr, z0, 0)
  xlim <- range(Re(allz)) + c(-1, 1) * 0.08 * diff(range(Re(allz)))
  ylim <- range(Im(allz)) + c(-1, 1) * 0.08 * diff(range(Im(allz)))
  graphics::plot(0, 0, type = "n", asp = 1, xlim = xlim, ylim = ylim,
                 xlab = "Real", ylab = "Imaginary",
                 main = sprintf("GSTATIS trajectories, view %d, dimension %d", view, dim), ...)
  graphics::abline(h = 0, v = 0, col = "grey80", lty = 3)
  nT <- ncol(tr)
  cols <- grDevices::hcl.colors(nrow(tr), "Dark 3")
  for (j in seq_len(nrow(tr))) {
    graphics::lines(Re(tr[j, ]), Im(tr[j, ]), col = cols[j])
    graphics::points(Re(tr[j, ]), Im(tr[j, ]), col = cols[j], pch = 1, cex = 0.6)
    graphics::text(Re(tr[j, nT]), Im(tr[j, nT]), labels = colnames(tr)[nT],
                   col = cols[j], cex = 0.5, pos = 4)
  }
  graphics::points(Re(z0), Im(z0), pch = 16, col = cols)
  graphics::text(Re(z0), Im(z0), labels = rownames(tr), pos = 3, cex = 0.8, col = cols)
  return(invisible(x))
}
