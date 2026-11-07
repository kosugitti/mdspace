# mdspace の共通の骨組み ------------------------------------------------------
# 全モデルの返り値は c("mdspace_{手法}", "mdspace") のクラスを持つ list。
# 共通スロット: conf, stress, call, niter, data, labels, model

# 返り値を組み立てる -----------------------------------------------------------
.new_mdspace <- function(model, conf, stress = NA_real_, niter = NA_integer_,
                         data = NULL, call = NULL, ...) {
  labels <- rownames(conf)
  if (is.null(labels)) {
    labels <- as.character(seq_len(nrow(conf)))
    rownames(conf) <- labels
  }
  obj <- list(
    model = model,
    conf = conf,
    stress = stress,
    niter = niter,
    data = data,
    labels = labels,
    call = call
  )
  obj <- c(obj, list(...))
  class(obj) <- c(paste0("mdspace_", model), "mdspace")
  return(obj)
}

# 共通の print -----------------------------------------------------------------
#' Print an mdspace object
#'
#' @param x An object of class `mdspace`.
#' @param digits Number of digits to print.
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.mdspace <- function(x, digits = 4, ...) {
  cat("mdspace:", x$model, "\n")
  if (!is.null(x$call)) {
    cat("Call: ", deparse(x$call), "\n", sep = "")
  }
  ndim <- if (is.null(x$ndim)) ncol(x$conf) else x$ndim
  cat("Objects:", nrow(x$conf), " Dimensions:", ndim, "\n")
  if (!is.na(x$stress)) {
    cat("Stress:", format(round(x$stress, digits), nsmall = digits), "\n")
  }
  if (!is.na(x$niter)) {
    cat("Iterations:", x$niter, "\n")
  }
  return(invisible(x))
}

# 共通の summary ---------------------------------------------------------------
#' Summarize an mdspace object
#'
#' @param object An object of class `mdspace`.
#' @param ... Not used.
#' @return An object of class `summary.mdspace`.
#' @export
summary.mdspace <- function(object, ...) {
  out <- list(
    model = object$model,
    call = object$call,
    stress = object$stress,
    niter = object$niter,
    conf = object$conf
  )
  class(out) <- "summary.mdspace"
  return(out)
}

#' @export
print.summary.mdspace <- function(x, digits = 4, ...) {
  cat("mdspace:", x$model, "\n")
  if (!is.na(x$stress)) {
    cat("Stress:", format(round(x$stress, digits), nsmall = digits), "\n")
  }
  cat("\nConfiguration:\n")
  print(round(x$conf, digits))
  return(invisible(x))
}

# 共通の plot: 布置の散布図 ------------------------------------------------------
#' Plot the configuration of an mdspace object
#'
#' @param x An object of class `mdspace`.
#' @param dims Two dimensions to plot.
#' @param labels Logical; draw object labels.
#' @param ... Passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
plot.mdspace <- function(x, dims = c(1, 2), labels = TRUE, ...) {
  conf <- x$conf
  if (ncol(conf) == 1) {
    conf <- cbind(conf, 0)
    dims <- c(1, 2)
  }
  xy <- conf[, dims, drop = FALSE]
  lim <- range(xy) * 1.15
  graphics::plot(xy, asp = 1, xlim = lim, ylim = lim,
                 xlab = paste("Dimension", dims[1]),
                 ylab = paste("Dimension", dims[2]),
                 pch = 16, main = x$model, ...)
  graphics::abline(h = 0, v = 0, col = "grey80", lty = 3)
  if (labels) {
    graphics::text(xy, labels = x$labels, pos = 3, cex = 0.8)
  }
  return(invisible(x))
}
