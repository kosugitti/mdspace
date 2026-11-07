# 近接性データの受け口と非対称の分解 ----------------------------------------------

#' Coerce data to a square proximity matrix
#'
#' Accepts a square matrix, a data frame, or a `dist` object and returns a
#' square numeric matrix with row and column labels. Rows are "from" objects
#' and columns are "to" objects, as in asymmetric MDS.
#'
#' @param x A square matrix, data frame, or `dist` object.
#' @param type Either `"dissimilarity"` or `"similarity"`. Stored as an
#'   attribute; the values are not transformed.
#' @param diag Value to put on the diagonal. `NULL` keeps the diagonal as is.
#' @return A square numeric matrix with attribute `type`.
#' @examples
#' d <- as_proximity(eurodist)
#' attr(d, "type")
#' @export
as_proximity <- function(x, type = c("dissimilarity", "similarity"),
                         diag = NULL) {
  type <- match.arg(type)
  if (inherits(x, "dist")) {
    x <- as.matrix(x)
  }
  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop("`x` must be numeric.")
  }
  if (nrow(x) != ncol(x)) {
    stop("`x` must be a square matrix.")
  }
  # ラベルを補う(行名と列名のどちらかがあれば両方にそろえる)
  lab <- rownames(x)
  if (is.null(lab)) lab <- colnames(x)
  if (is.null(lab)) lab <- as.character(seq_len(nrow(x)))
  dimnames(x) <- list(lab, lab)
  if (!is.null(diag)) {
    base::diag(x) <- diag
  }
  attr(x, "type") <- type
  return(x)
}

#' Decompose a square matrix into symmetric and skew-symmetric parts
#'
#' `X = S + K` with `S = (X + t(X)) / 2` and `K = (X - t(X)) / 2`. The share of
#' asymmetry is `sum(K^2) / sum(X^2)` (Gower, 1977). Missing values are ignored
#' in the sums.
#'
#' @param x A square proximity matrix (see [as_proximity()]).
#' @return An object of class `mdspace_asym` with elements `S`, `K`, `ss`
#'   (sums of squares of `X`, `S`, `K`) and `asymmetry`.
#' @references Gower, J. C. (1977). The analysis of asymmetry and
#'   orthogonality. In J. R. Barra et al. (Eds.), *Recent developments in
#'   statistics* (pp. 109-123). North-Holland.
#' @examples
#' x <- matrix(c(0, 1, 1, 7, 1, 0, 1, 7, 7, 7, 0, 1, 1, 1, 7, 0), 4, byrow = TRUE)
#' asym_decomp(x)
#' @export
asym_decomp <- function(x) {
  x <- as_proximity(x)
  S <- (x + t(x)) / 2
  K <- (x - t(x)) / 2
  attr(S, "type") <- NULL
  attr(K, "type") <- NULL
  ss <- c(total = sum(x^2, na.rm = TRUE),
          symmetric = sum(S^2, na.rm = TRUE),
          skew = sum(K^2, na.rm = TRUE))
  out <- list(S = S, K = K, ss = ss,
              asymmetry = unname(ss["skew"] / ss["total"]))
  class(out) <- "mdspace_asym"
  return(out)
}

#' @export
print.mdspace_asym <- function(x, digits = 4, ...) {
  cat("Symmetric / skew-symmetric decomposition\n")
  cat("Sum of squares: total =", format(x$ss["total"], digits = digits),
      " symmetric =", format(x$ss["symmetric"], digits = digits),
      " skew =", format(x$ss["skew"], digits = digits), "\n")
  cat("Share of asymmetry:", format(round(x$asymmetry, digits), nsmall = digits), "\n")
  return(invisible(x))
}

#' Plot the skew-symmetric part as a heat map
#'
#' Positive cells (row > column) are red and negative cells are blue.
#'
#' @param x An object of class `mdspace_asym`.
#' @param ... Passed to [graphics::image()].
#' @return `x`, invisibly.
#' @export
plot.mdspace_asym <- function(x, ...) {
  K <- x$K
  n <- nrow(K)
  m <- max(abs(K), na.rm = TRUE)
  if (m == 0) m <- 1
  cols <- grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(21)
  # image() は行列を転置して下から描くので，表と同じ向きになるよう並べ替える
  graphics::image(seq_len(n), seq_len(n), t(K[n:1, , drop = FALSE]),
                  col = cols, zlim = c(-m, m), axes = FALSE,
                  xlab = "to", ylab = "from", main = "Skew-symmetric part", ...)
  graphics::axis(1, at = seq_len(n), labels = colnames(K), las = 2, cex.axis = 0.8)
  graphics::axis(2, at = seq_len(n), labels = rev(rownames(K)), las = 2, cex.axis = 0.8)
  graphics::box()
  return(invisible(x))
}
