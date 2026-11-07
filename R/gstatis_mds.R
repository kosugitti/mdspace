# GSTATIS(二相三元の HFM; Chino, Grorud & Yoshino, 1996) --------------------------

#' GSTATIS: STATIS for a set of asymmetric matrices via Hermitian forms
#'
#' A set of `T` square asymmetric matrices on the same `n` objects (e.g. the
#' same group rated at `T` occasions) is analyzed in three steps.
#'
#' 1. **Interstructure.** Each matrix is turned into its Hermitian form
#'    `H_t = S_t + iK_t`. The generalized RV coefficient between occasions is
#'    `Re tr(H_s H_t) / sqrt(tr(H_s H_s) tr(H_t H_t))`, and the `T x T` RV
#'    matrix is eigen-decomposed.
#' 2. **Compromise.** The `k`-th eigenvector `v_k` of the RV matrix gives the
#'    weights of a compromise `C_k = sum_t v_tk H_t`. The first one is the
#'    usual compromise; later ones are other "points of view". Each compromise
#'    is analyzed by HFM (see [hfm_mds()]).
#' 3. **Trajectories.** Each occasion is projected onto the planes of the first
#'    compromise: `Z_t = H_t U |Lambda|^{-1/2} sign(Lambda)`, where `U` and
#'    `Lambda` are the eigenvectors and eigenvalues of `C_1`. If `H_t` equals
#'    `C_1`, `Z_t` equals the compromise coordinates.
#'
#' The first eigenvector of the RV matrix is signed so that its weights sum to
#' a positive value and then rescaled to sum to 1, so that the first
#' compromise is a weighted mean of the occasions and the trajectories are on
#' the same scale as the compromise. (The 2001 program used the raw
#' eigenvector; the configuration differs only by a constant factor.) Later
#' points of view use the raw eigenvectors. With `scale = "norm"` each `H_t` is divided by its
#' Frobenius norm before step 2 (the usual STATIS normalization); the default
#' `"none"` follows the 2001 Visual Basic program and the 2014 R code.
#'
#' @param x A list of square matrices with the same labels, or an
#'   `n x n x T` array.
#' @param ndim Number of complex dimensions to keep.
#' @param ncomp Number of compromises (points of view) to build. Use
#'   `ncomp = NULL` with `threshold` to keep every RV component whose
#'   eigenvalue exceeds `threshold`, as in the 2001 program.
#' @param threshold Eigenvalue threshold used when `ncomp = NULL`.
#' @param scale `"none"` or `"norm"`; see Details.
#' @param rotate Rotation of each complex plane; see [hfm_mds()].
#' @return An object of class `c("mdspace_gstatis", "mdspace")` with
#'   `conf`, `coord`, `eigen`, `gof`, `unexplained` (as in [hfm_mds()], for the
#'   first compromise), `rv` (RV matrix), `rv_eigen` (its eigenvalues and
#'   vectors), `weights` (`T x ncomp`), `compromise` (list of HFM fits, one per
#'   point of view), and `trajectory` (complex array `n x ndim x T`).
#' @references
#' Chino, N., Grorud, A., & Yoshino, R. (1996). A complex analysis for
#'   two-mode three-way asymmetric relational data. Paper presented at the
#'   Fifth Conference of the International Federation of Classification
#'   Societies, Kobe, Japan.
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
                        rotate = c("std", "varimax", "none")) {
  cl <- match.call()
  scale <- match.arg(scale)
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
  if (any(vapply(X, anyNA, logical(1)))) {
    stop("gstatis_mds() does not allow missing values.")
  }
  labs <- rownames(X[[1]])
  X <- lapply(X, function(m) {
    dimnames(m) <- list(labs, labs)
    return(m)
  })
  occ <- names(x)
  if (is.null(occ)) occ <- paste0("T", seq_along(X))
  nT <- length(X)
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
  if (sum(V[, 1]) < 0) {
    V[, 1] <- -V[, 1]
  }
  # 第1の妥協は重みの和を 1 にして時点の加重平均にする
  W <- V
  W[, 1] <- V[, 1] / sum(V[, 1])
  dimnames(W) <- list(occ, paste0("C", seq_len(nT)))
  dimnames(V) <- list(occ, paste0("C", seq_len(nT)))

  if (is.null(ncomp)) {
    ncomp <- max(1, sum(rv_eig$values > threshold))
  }
  ncomp <- min(ncomp, nT)

  # 第2段: 視点ごとの妥協行列を HFM にかける ------------------------------------------
  comp <- vector("list", ncomp)
  for (k in seq_len(ncomp)) {
    Ck <- Reduce(`+`, Map(function(h, w) h * w, H, W[, k]))
    eg <- .hfm_eigen(Ck, rotate = rotate)
    lam <- eg$values
    keep <- seq_len(ndim)
    coord <- eg$vectors[, keep, drop = FALSE] %*% diag(sqrt(abs(lam[keep])), ndim)
    dimnames(coord) <- list(labs, paste0("D", keep))
    gof <- lam^2 / sum(lam^2)
    comp[[k]] <- list(matrix = Ck, vectors = eg$vectors, eigen = lam,
                      coord = coord, gof = gof,
                      unexplained = 1 - sum(gof[keep]))
  }
  names(comp) <- paste0("C", seq_len(ncomp))

  # 第3段: 各時点を第1の妥協空間へ射影して軌跡を作る -----------------------------------
  c1 <- comp[[1]]
  keep <- seq_len(ndim)
  U1 <- c1$vectors[, keep, drop = FALSE]
  lam1 <- c1$eigen[keep]
  if (any(abs(lam1) < 1e-10 * max(abs(c1$eigen)))) {
    warning("Some kept eigenvalues of the compromise are near 0; trajectories are unstable.")
  }
  P <- U1 %*% diag(sign(lam1) / sqrt(abs(lam1)), ndim)
  traj <- array(0i, dim = c(n, ndim, nT),
                dimnames = list(labs, paste0("D", keep), occ))
  for (t in seq_len(nT)) {
    traj[, , t] <- H[[t]] %*% P
  }

  out <- .new_mdspace("gstatis", conf = .complex_to_conf(c1$coord), data = X,
                      call = cl, ndim = ndim, coord = c1$coord,
                      eigen = c1$eigen, gof = c1$gof,
                      unexplained = c1$unexplained,
                      rv = rv, rv_eigen = list(values = rv_eig$values, vectors = V),
                      weights = W[, seq_len(ncomp), drop = FALSE],
                      compromise = comp, trajectory = traj, scale = scale)
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
  cat("\nWeights of the first compromise:\n")
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
#' `type = "interstructure"` plots the occasions on the first two eigenvectors
#' of the RV matrix. `type = "compromise"` plots one complex plane of a
#' compromise like [plot.mdspace_hfm()]. `type = "trajectory"` draws, for each
#' object, its positions over the occasions joined by lines, with the
#' compromise position as a filled point.
#'
#' @param x An object of class `mdspace_gstatis`.
#' @param type `"trajectory"`, `"compromise"`, or `"interstructure"`.
#' @param dim Which complex dimension to plot.
#' @param view Which compromise (point of view) for `type = "compromise"`.
#' @param ... Passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
plot.mdspace_gstatis <- function(x, type = c("trajectory", "compromise", "interstructure"),
                                 dim = 1, view = 1, ...) {
  type <- match.arg(type)
  if (type == "interstructure") {
    V <- x$rv_eigen$vectors
    ev <- x$rv_eigen$values
    xy <- V[, 1:2, drop = FALSE] %*% diag(sqrt(pmax(ev[1:2], 0)), 2)
    lim <- c(-1, 1) * max(abs(xy)) * 1.2
    graphics::plot(xy, asp = 1, xlim = lim, ylim = lim, pch = 16,
                   xlab = "Component 1", ylab = "Component 2",
                   main = "GSTATIS interstructure", ...)
    graphics::abline(h = 0, v = 0, col = "grey80", lty = 3)
    graphics::text(xy, labels = rownames(V), pos = 3, cex = 0.8)
    return(invisible(x))
  }
  if (type == "compromise") {
    cp <- x$compromise[[view]]
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
  # 軌跡: 時点ごとの位置を線でつなぎ，妥協の位置を黒丸で示す
  tr <- x$trajectory[, dim, , drop = TRUE]
  z0 <- x$coord[, dim]
  # 原点も含めた点の範囲に合わせて描く(縦横比は 1)
  allz <- c(tr, z0, 0)
  xlim <- range(Re(allz)) + c(-1, 1) * 0.08 * diff(range(Re(allz)))
  ylim <- range(Im(allz)) + c(-1, 1) * 0.08 * diff(range(Im(allz)))
  graphics::plot(0, 0, type = "n", asp = 1, xlim = xlim, ylim = ylim,
                 xlab = "Real", ylab = "Imaginary",
                 main = sprintf("GSTATIS trajectories, dimension %d", dim), ...)
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
  graphics::text(Re(z0), Im(z0), labels = x$labels, pos = 3, cex = 0.8, col = cols)
  return(invisible(x))
}
