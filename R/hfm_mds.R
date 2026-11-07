# HFM(Hermitian Form Model; Chino & Shiraiwa, 1993) ----------------------------

#' Hermitian form model for asymmetric proximities
#'
#' The square matrix `X` is split into its symmetric part `S` and its
#' skew-symmetric part `K`, and the Hermitian matrix `H = S + iK` is
#' eigen-decomposed: `H = U diag(lambda) U*`. The eigenvalues are real and are
#' ordered by absolute value. Each dimension is a complex plane; object `j` is
#' the point `z_j = u_j * sqrt(|lambda|)`. For a dimension with a positive
#' eigenvalue, `Re(z_j Conj(z_k))` reproduces the symmetric part and
#' `Im(z_j Conj(z_k))` the skew-symmetric part, so the angle from `k` to `j`
#' (counterclockwise) shows the direction of asymmetry.
#'
#' An eigenvector is determined only up to a factor `exp(i * theta)`, i.e. a
#' rotation of its complex plane. `rotate = "std"` (default) turns each plane
#' so that the principal axis of the points lies on the real axis (the sum of
#' squared real parts is maximal). `rotate = "varimax"` turns each plane so
#' that the variance of `Re(z^2)` over objects is maximal, which pushes the
#' points toward the real and imaginary axes. Both follow the 2001 Visual Basic
#' program (`rotateSTD` and `rotateMAX`). In both cases the plane is then
#' flipped, if needed, so that the object with the largest absolute real part
#' lies on the positive side.
#'
#' @param x A square similarity matrix (rows = from, columns = to). The
#'   diagonal is used as is.
#' @param ndim Number of complex dimensions to keep.
#' @param rotate `"std"`, `"varimax"`, or `"none"`.
#' @return An object of class `c("mdspace_hfm", "mdspace")` with
#'   `coord` (n x ndim complex matrix), `conf` (n x 2*ndim real matrix: real and
#'   imaginary parts of each dimension), `eigen` (all eigenvalues), `gof`
#'   (`lambda^2 / sum(lambda^2)` for each dimension), `unexplained` (share of
#'   the sum of squares of `H` not reproduced by the kept dimensions, i.e.
#'   1 - cumulative GOF), and `hermitian` (`H`). HFM is not a distance model,
#'   so `stress` is `NA`.
#' @references Chino, N., & Shiraiwa, K. (1993). Geometrical structures of
#'   some non-distance models for asymmetric MDS. *Behaviormetrika*, 20(1),
#'   35-47.
#' @examples
#' x <- matrix(c(0, 1, 1, 7, 1, 0, 1, 7, 7, 7, 0, 1, 1, 1, 7, 0), 4, byrow = TRUE)
#' fit <- hfm_mds(x, ndim = 2)
#' fit
#' plot(fit, dim = 1)
#' @export
hfm_mds <- function(x, ndim = 2, rotate = c("std", "varimax", "none")) {
  cl <- match.call()
  rotate <- match.arg(rotate)
  X <- as_proximity(x, type = "similarity")
  if (anyNA(X)) {
    stop("hfm_mds() does not allow missing values.")
  }
  n <- nrow(X)
  if (ndim < 1 || ndim > n) {
    stop("`ndim` must be between 1 and n.")
  }

  # エルミート行列 H = S + iK を固有値分解し，各平面の回転を固定する ------------------
  H <- .hermitian(X)
  eg <- .hfm_eigen(H, rotate = rotate)
  lam <- eg$values
  keep <- seq_len(ndim)
  coord <- eg$vectors[, keep, drop = FALSE] %*% diag(sqrt(abs(lam[keep])), ndim)
  rownames(coord) <- rownames(X)
  colnames(coord) <- paste0("D", keep)
  conf <- .complex_to_conf(coord)

  # 適合の指標: 次元ごとの GOF と，残した次元で再現できなかった二乗和の割合
  # (距離モデルではないので stress は使わない)
  gof <- lam^2 / sum(lam^2)
  unexplained <- 1 - sum(gof[keep])

  out <- .new_mdspace("hfm", conf = conf, data = X,
                      call = cl, ndim = ndim, coord = coord, eigen = lam, gof = gof, unexplained = unexplained,
                      hermitian = H)
  return(out)
}

#' @export
print.mdspace_hfm <- function(x, digits = 4, ...) {
  NextMethod()
  k <- ncol(x$coord)
  cat("Unexplained SS (1 - cumulative GOF):",
      format(round(x$unexplained, digits), nsmall = digits), "\n")
  tab <- rbind(eigenvalue = x$eigen[seq_len(k)], gof = x$gof[seq_len(k)],
               cumulative = cumsum(x$gof[seq_len(k)]))
  colnames(tab) <- colnames(x$coord)
  cat("\n")
  print(round(tab, digits))
  if (any(x$eigen[seq_len(k)] < 0)) {
    cat("Note: negative eigenvalues reverse the reading of angles in those planes.\n")
  }
  return(invisible(x))
}

#' Plot one complex plane of an HFM solution
#'
#' Each object is drawn as an arrow from the origin. For a positive
#' eigenvalue, the counterclockwise angle from object k to object j is
#' proportional to the skew-symmetric part `K[j, k]`.
#'
#' @param x An object of class `mdspace_hfm`.
#' @param dim Which complex dimension to plot.
#' @param ... Passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
plot.mdspace_hfm <- function(x, dim = 1, ...) {
  z <- x$coord[, dim]
  lim <- c(-1, 1) * max(Mod(z)) * 1.2
  graphics::plot(Re(z), Im(z), asp = 1, xlim = lim, ylim = lim, type = "n",
                 xlab = "Real", ylab = "Imaginary",
                 main = sprintf("HFM dimension %d (eigenvalue %.3f)",
                                dim, x$eigen[dim]), ...)
  graphics::abline(h = 0, v = 0, col = "grey80", lty = 3)
  # 原点のごく近くの対象は矢印を描かない(短すぎる矢印は警告になる)
  ok <- Mod(z) > 1e-3 * max(Mod(z))
  graphics::arrows(0, 0, Re(z[ok]), Im(z[ok]), length = 0.08, col = "grey40")
  graphics::text(Re(z), Im(z), labels = x$labels, pos = 3, cex = 0.8)
  return(invisible(x))
}

# 1 本の複素ベクトルの回転を決める ------------------------------------------------
# std: sum(z^2) の偏角の半分だけ戻す＝主軸を実軸へ(Re の二乗和が最大)
# varimax: Re(z^2) の分散が最大になる角＝arg(sum(z^4) - sum(z^2)^2 / n) / 4
#   (90 度ずつの不定が残るので，その中で Re の二乗和が最大のものを選ぶ)
.rotate_plane <- function(z, method = c("std", "varimax")) {
  method <- match.arg(method)
  n <- length(z)
  if (method == "std") {
    th <- Arg(sum(z^2)) / 2
    cand <- th
  } else {
    w <- sum(z^4) - sum(z^2)^2 / n
    th <- Arg(w) / 4
    cand <- th + (0:3) * pi / 2
  }
  ss <- vapply(cand, function(a) sum(Re(z * exp(-1i * a))^2), numeric(1))
  z <- z * exp(-1i * cand[which.max(ss)])
  # 実部の絶対値が最大の対象を正の側に置く
  j <- which.max(abs(Re(z)))
  if (Re(z[j]) < 0) {
    z <- -z
  }
  return(z)
}

# 正方行列からエルミート行列 H = S + iK を作る --------------------------------------
.hermitian <- function(X) {
  dec <- asym_decomp(X)
  n <- nrow(X)
  H <- matrix(complex(real = dec$S, imaginary = dec$K), n, n,
              dimnames = dimnames(dec$S))
  return(H)
}

# エルミート行列の固有値分解: 絶対値の大きい順に並べ，各平面の回転を固定する ---------------
.hfm_eigen <- function(H, rotate = c("std", "varimax", "none")) {
  rotate <- match.arg(rotate)
  ev <- eigen(H, symmetric = TRUE)
  lam <- Re(ev$values)
  ord <- order(abs(lam), decreasing = TRUE)
  lam <- lam[ord]
  U <- ev$vectors[, ord, drop = FALSE]
  if (rotate != "none") {
    for (t in seq_len(ncol(U))) {
      U[, t] <- .rotate_plane(U[, t], method = rotate)
    }
  }
  return(list(values = lam, vectors = U))
}

# 複素座標(n x ndim)を実部・虚部を並べた実行列(n x 2ndim)にする -----------------------
.complex_to_conf <- function(coord) {
  ndim <- ncol(coord)
  conf <- matrix(0, nrow(coord), 2 * ndim)
  conf[, seq(1, 2 * ndim, 2)] <- Re(coord)
  conf[, seq(2, 2 * ndim, 2)] <- Im(coord)
  rownames(conf) <- rownames(coord)
  colnames(conf) <- paste0(rep(c("Re", "Im"), ndim), rep(seq_len(ndim), each = 2))
  return(conf)
}
