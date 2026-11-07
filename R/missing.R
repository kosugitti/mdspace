# 非対称行列の欠測: 型の集計と補完(MINERS 法) ---------------------------------------

#' Missing-value patterns of asymmetric proximity matrices
#'
#' Classifies the missing cells of a square matrix (or each matrix of a list)
#' into four types:
#'
#' * `one_side`: off-diagonal cell `(j, k)` is missing but `(k, j)` is
#'   observed (this includes an object that did not respond, whose row is
#'   missing while its column is observed);
#' * `both_sides`: both `(j, k)` and `(k, j)` are missing (counted once per
#'   cell, so a pair counts 2);
#' * `diagonal`: missing diagonal cells;
#' * `absent`: objects whose whole row and column (off the diagonal) are
#'   missing, i.e. objects absent at that occasion. Their cells are not
#'   counted in the other types.
#'
#' @param x A square matrix, or a list of square matrices.
#' @return A data frame with one row per matrix: `n`, `one_side`,
#'   `both_sides`, `diagonal`, `absent` (number of absent objects),
#'   `absent_objects` (their labels), and `rate` (share of missing cells).
#' @examples
#' x <- matrix(c(NA, 1, 2, NA, NA, 3, 4, 5, 6), 3, byrow = TRUE)
#' missing_pattern(x)
#' @export
missing_pattern <- function(x) {
  if (is.matrix(x) || is.data.frame(x)) {
    x <- list(x)
  }
  nm <- names(x)
  if (is.null(nm)) nm <- paste0("T", seq_along(x))
  rows <- lapply(seq_along(x), function(i) {
    p <- .missing_types(as_proximity(x[[i]]))
    return(data.frame(matrix = nm[i], n = p$n, one_side = p$one_side,
                      both_sides = p$both_sides, diagonal = p$diagonal,
                      absent = length(p$absent),
                      absent_objects = paste(p$absent, collapse = ","),
                      rate = p$rate))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  return(out)
}

# 欠測の型を数える(内部) ------------------------------------------------------------
.missing_types <- function(X) {
  n <- nrow(X)
  m <- is.na(X)
  off <- row(m) != col(m)
  # その時点にいない対象: 対角以外の行と列がすべて欠測
  absent_idx <- which(rowSums(m & off) == n - 1 & colSums(m & off) == n - 1)
  in_absent <- matrix(FALSE, n, n)
  in_absent[absent_idx, ] <- TRUE
  in_absent[, absent_idx] <- TRUE
  rest <- m & !in_absent
  return(list(
    n = n,
    one_side = sum(rest & off & !t(m)),
    both_sides = sum(rest & off & t(m)),
    diagonal = sum(diag(rest)),
    absent = rownames(X)[absent_idx],
    rate = mean(m)
  ))
}

# MINERS 法による補完(内部) ---------------------------------------------------------
# 1. 片側だけの欠けは転置の値で埋めて固定する(その組の歪対称部は 0 とみなす)
# 2. 対角と両側の欠けは，行と列の観測値の平均を初期値にして，
#    HFM の上位 p 次元(絶対値の大きい順)で再構成した値に置き換えることを収束まで繰り返す
#    (再構成した H から X = Re(H) + Im(H) で元の行列の値に戻す)
.impute_miners <- function(X, p, tol = 1e-4, maxit = 1000, label = NULL) {
  types <- .missing_types(X)
  if (length(types$absent) > 0) {
    stop(sprintf("Object(s) %s %s have no observed proximities (absent); remove them or the occasion.",
                 paste(types$absent, collapse = ", "),
                 if (is.null(label)) "" else paste0("at ", label)))
  }
  m <- is.na(X)
  W <- X
  off <- row(m) != col(m)
  one <- m & off & !t(m)
  W[one] <- t(X)[one]
  iter_cells <- m & !one
  type <- ifelse(row(m) == col(m), "diagonal", ifelse(one, "one_side", "both_sides"))
  niter <- 0L
  converged <- TRUE
  if (any(iter_cells)) {
    # 初期値: 行 j と列 k の観測値の平均
    for (k in which(iter_cells)) {
      j <- row(m)[k]
      l <- col(m)[k]
      W[k] <- mean(c(X[j, ], X[, l]), na.rm = TRUE)
    }
    converged <- FALSE
    for (it in seq_len(maxit)) {
      eg <- .hfm_eigen(.hermitian(W), rotate = "none")
      keep <- seq_len(p)
      U <- eg$vectors[, keep, drop = FALSE]
      Hh <- U %*% diag(eg$values[keep], p) %*% Conj(t(U))
      R <- Re(Hh) + Im(Hh)
      change <- max(abs(R[iter_cells] - W[iter_cells]))
      W[iter_cells] <- R[iter_cells]
      niter <- it
      if (change < tol) {
        converged <- TRUE
        break
      }
    }
    if (!converged) {
      warning(sprintf("MINERS imputation did not converge in %d iterations%s.", maxit,
                      if (is.null(label)) "" else paste0(" (", label, ")")))
    }
  }
  idx <- which(m, arr.ind = TRUE)
  imputed <- data.frame(row = rownames(X)[idx[, 1]], col = colnames(X)[idx[, 2]],
                        type = type[idx], value = W[idx])
  attr(W, "type") <- attr(X, "type")
  return(list(data = W, imputed = imputed, niter = niter, converged = converged,
              pattern = types))
}
