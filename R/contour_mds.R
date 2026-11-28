# 等高線マップモデル(contour map model): Abelson (1954-55)・小杉・藤原 (2004) ----------

#' Contour map model of valences on an MDS configuration
#'
#' Places evaluations of the objects (for example, liking ratings) on an MDS
#' configuration as a field of force (Abelson, 1954-55; Kosugi & Fujihara,
#' 2004). With `v_j` the centered rating of object `j` and `d(P, x_j)` the
#' Euclidean distance from a point `P` to object `j` in the configuration, the
#' valence at `P` is
#' `V(P) = sum_j v_j / (1 + d(P, x_j)^2)`.
#' It is the potential of point charges `v_j` with `d` replaced by
#' `1 + d^2`, so that it stays finite on the objects. The contour line of
#' valence 0 divides the space into liked and disliked regions; steep slopes
#' and several hills are read as tension and as a lack of integration of the
#' attitude (Abelson, 1954-55).
#'
#' The height of the surface at object `j` is
#' `V(x_j) = v_j + sum_{k != j} v_k / (1 + d_jk^2)`, so the discrepancy
#' `v_j - V(x_j)` (the "length of the perpendicular" in Kosugi & Fujihara,
#' 2004) is minus the force that the other objects exert on `j`. It is large
#' for an object rated in the opposite direction to the objects near it.
#'
#' The constant 1 in the denominator fixes the unit of length, so the map
#' depends on the scale of the configuration: enlarging it makes the hills
#' sharper. By default the configuration is used as given, as in the
#' original. `normalize = TRUE` first centers it and scales it so that the
#' mean squared distance of the objects from their centroid is 1.
#'
#' Ratings that are `NA` are left out of the sum (their charge is 0) and
#' their discrepancy is `NA`.
#'
#' @param conf A configuration (n x p matrix, p >= 2) or an `mdspace` object,
#'   whose `$conf` is used.
#' @param ratings A numeric vector of length n (one person) or a matrix with
#'   one row per person and one column per object. When both `ratings` and
#'   the configuration have object names, they are matched by name.
#' @param center Value subtracted from the ratings so that 0 is neutral (for
#'   example 3 for a five-point scale, as in Kosugi & Fujihara, 2004).
#' @param normalize Logical; scale the configuration so that the mean squared
#'   distance of the objects from their centroid is 1.
#' @return An object of class `c("mdspace_contour", "mdspace")` with
#'   \describe{
#'     \item{conf}{The configuration used (after `normalize`).}
#'     \item{ratings}{The centered ratings (persons x objects).}
#'     \item{height}{The valence at each object (persons x objects).}
#'     \item{discrepancy}{`ratings - height` (persons x objects).}
#'     \item{center, normalize}{The arguments.}
#'   }
#'   `stress` and `niter` are `NA` because nothing is estimated.
#' @references
#' Abelson, R. P. (1954-55). A technique and a model for multi-dimensional
#' scaling. *Public Opinion Quarterly*, *18*(4), 405-418.
#'
#' Kosugi, K., & Fujihara, T. (2004). A technique for mapping attitude by a
#' contour map model. *Behaviormetrika (Kodo Keiryogaku)*, *31*(1), 17-24.
#' (in Japanese)
#' @seealso [predict.mdspace_contour()] for valences at any point,
#'   [plot.mdspace_contour()] for the contour map and the valence surface.
#' @examples
#' conf <- matrix(c(2.04, -0.48, -0.30, 0.43, 1.08, 1.04,
#'                  -1.66, 0.54, -0.41, -0.53, -0.75, -1.01),
#'                ncol = 2, byrow = TRUE,
#'                dimnames = list(c("Giants", "Dragons", "Tigers", "Carp",
#'                                  "Swallows", "Baystars"), NULL))
#' fit <- contour_mds(conf, c(1, 2, 1, 5, 2, 2), center = 3)
#' summary(fit)
#' plot(fit)
#' plot(fit, type = "persp")
#' @export
contour_mds <- function(conf, ratings, center = 0, normalize = FALSE) {
  cl <- match.call()
  if (inherits(conf, "mdspace")) {
    conf <- conf$conf
  }
  conf <- as.matrix(conf)
  if (!is.numeric(conf) || anyNA(conf)) {
    stop("`conf` must be a numeric matrix without missing values.")
  }
  n <- nrow(conf)
  if (ncol(conf) < 2) {
    stop("`conf` must have at least 2 dimensions.")
  }
  if (is.null(rownames(conf))) {
    rownames(conf) <- as.character(seq_len(n))
  }
  if (!is.numeric(center) || length(center) != 1 || is.na(center)) {
    stop("`center` must be a single number.")
  }
  if (normalize) {
    conf <- sweep(conf, 2, colMeans(conf))
    conf <- conf * sqrt(n / sum(conf^2))
  }
  R <- .contour_ratings(ratings, rownames(conf)) - center
  # 対象の点での高さ(自分の項 v_j も含む)
  W <- 1 / (1 + as.matrix(stats::dist(conf))^2)
  R0 <- R
  R0[is.na(R0)] <- 0
  H <- R0 %*% W
  dimnames(H) <- dimnames(R)
  out <- .new_mdspace(model = "contour", conf = conf, data = ratings, call = cl,
                      ratings = R, height = H, discrepancy = R - H,
                      center = center, normalize = normalize)
  return(out)
}

# 評定を人 x 対象の行列にそろえる(内部) ------------------------------------------
.contour_ratings <- function(ratings, labels) {
  n <- length(labels)
  if (is.data.frame(ratings)) {
    ratings <- as.matrix(ratings)
  }
  if (is.null(dim(ratings))) {
    nm <- names(ratings)
    ratings <- matrix(ratings, nrow = 1, dimnames = list("1", nm))
  }
  if (!is.numeric(ratings)) {
    stop("`ratings` must be numeric.")
  }
  if (!is.null(colnames(ratings)) && !all(labels == as.character(seq_len(n)))) {
    miss <- setdiff(labels, colnames(ratings))
    if (length(miss) > 0 || ncol(ratings) != n) {
      stop(sprintf("The names of `ratings` do not match the objects of `conf` (missing: %s).",
                   paste(miss, collapse = ", ")))
    }
    ratings <- ratings[, labels, drop = FALSE]
  } else if (ncol(ratings) != n) {
    stop(sprintf("`ratings` must have %d values (one per object).", n))
  }
  colnames(ratings) <- labels
  if (is.null(rownames(ratings))) {
    rownames(ratings) <- as.character(seq_len(nrow(ratings)))
  }
  if (any(rowSums(!is.na(ratings)) == 0)) {
    stop("Each person must have at least one observed rating.")
  }
  return(ratings)
}

# 任意の点の誘因価(内部): P は m x p，返り値は m x 人 ---------------------------
.contour_valence <- function(P, conf, R) {
  D2 <- 0
  for (k in seq_len(ncol(conf))) {
    D2 <- D2 + outer(P[, k], conf[, k], "-")^2
  }
  R0 <- R
  R0[is.na(R0)] <- 0
  return((1 / (1 + D2)) %*% t(R0))
}

# 人を番号か名前で選ぶ(内部) ----------------------------------------------------
.contour_person <- function(object, person) {
  ids <- rownames(object$ratings)
  if (is.character(person)) {
    k <- match(person, ids)
  } else {
    k <- person
  }
  if (length(k) != 1 || is.na(k) || k < 1 || k > length(ids)) {
    stop("`person` must be one row number or name of the ratings.")
  }
  return(k)
}

#' Valences of a contour map model at given points
#'
#' @param object An object of class `mdspace_contour`.
#' @param newdata A matrix of points with as many columns as the
#'   configuration. `NULL` returns the heights at the objects.
#' @param ... Not used.
#' @return A matrix with one row per point and one column per person.
#' @export
predict.mdspace_contour <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    return(t(object$height))
  }
  P <- as.matrix(newdata)
  if (is.null(dim(newdata))) {
    P <- matrix(newdata, nrow = 1)
  }
  if (ncol(P) != ncol(object$conf)) {
    stop(sprintf("`newdata` must have %d columns.", ncol(object$conf)))
  }
  V <- .contour_valence(P, object$conf, object$ratings)
  colnames(V) <- rownames(object$ratings)
  return(V)
}

#' @export
print.mdspace_contour <- function(x, digits = 4, ...) {
  NextMethod()
  cat("Persons:", nrow(x$ratings), " Center:", x$center,
      " Normalized configuration:", x$normalize, "\n")
  return(invisible(x))
}

#' @export
summary.mdspace_contour <- function(object, ...) {
  out <- NextMethod()
  if (nrow(object$ratings) == 1) {
    out$extra <- list(
      "Rating, height and discrepancy" =
        data.frame(object = colnames(object$ratings), rating = object$ratings[1, ],
                   height = object$height[1, ], discrepancy = object$discrepancy[1, ])
    )
  } else {
    out$extra <- list("Height (persons x objects)" = object$height,
                      "Discrepancy (persons x objects)" = object$discrepancy)
  }
  return(out)
}

#' Plot a contour map model
#'
#' Draws the valences of one person over a plane of the configuration, as a
#' contour map (`type = "contour"`) with the contour of valence 0 drawn
#' thick, or as a valence surface (`type = "persp"`) with each object drawn at
#' its rating, filled when the rating is above the surface and open when it
#' is below (as in Kosugi & Fujihara, 2004). For three or more dimensions the
#' plane passes through the centroid of the configuration.
#'
#' @param x An object of class `mdspace_contour`.
#' @param person Row number or name of the person to draw.
#' @param type `"contour"` or `"persp"`.
#' @param dims Two dimensions spanning the plane.
#' @param grid Number of grid points on each axis. `NULL` uses 101 for
#'   `type = "contour"` and 41 for `type = "persp"`.
#' @param lim Range of both axes. `NULL` extends the range of the
#'   configuration by 30%.
#' @param nlevels Approximate number of contour levels.
#' @param labels Logical; draw object labels.
#' @param theta,phi Viewing angles for `type = "persp"`.
#' @param ... Passed to [graphics::contour()] or [graphics::persp()].
#' @return The grid (`x`, `y`, `z`), invisibly.
#' @export
plot.mdspace_contour <- function(x, person = 1, type = c("contour", "persp"),
                                 dims = c(1, 2), grid = NULL, lim = NULL, nlevels = 10,
                                 labels = TRUE, theta = 30, phi = 30, ...) {
  type <- match.arg(type)
  if (is.null(grid)) {
    grid <- if (type == "contour") 101 else 41
  }
  k <- .contour_person(x, person)
  conf <- x$conf
  if (length(dims) != 2 || any(!dims %in% seq_len(ncol(conf)))) {
    stop("`dims` must be two dimensions of the configuration.")
  }
  if (is.null(lim)) {
    r <- range(conf[, dims])
    lim <- r + c(-1, 1) * 0.3 * diff(r)
  }
  gx <- seq(lim[1], lim[2], length.out = grid)
  P <- matrix(colMeans(conf), nrow = grid * grid, ncol = ncol(conf), byrow = TRUE)
  P[, dims[1]] <- rep(gx, times = grid)
  P[, dims[2]] <- rep(gx, each = grid)
  z <- matrix(.contour_valence(P, conf, x$ratings[k, , drop = FALSE]), grid, grid)
  main <- if (nrow(x$ratings) == 1) "Contour map model" else
    paste("Contour map model:", rownames(x$ratings)[k])
  xlab <- paste("Dimension", dims[1])
  ylab <- paste("Dimension", dims[2])
  xy <- conf[, dims, drop = FALSE]
  if (type == "contour") {
    lev <- pretty(range(z), nlevels)
    graphics::contour(gx, gx, z, levels = setdiff(lev, 0), asp = 1, xlab = xlab, ylab = ylab,
                      main = main, col = "grey40", ...)
    if (min(z) < 0 && max(z) > 0) {
      graphics::contour(gx, gx, z, levels = 0, add = TRUE, lwd = 2, col = "firebrick")
    }
    graphics::points(xy, pch = 16)
    if (labels) {
      graphics::text(xy, labels = x$labels, pos = 3, cex = 0.8)
    }
  } else {
    v <- x$ratings[k, ]
    zl <- range(z, v, na.rm = TRUE)
    pm <- graphics::persp(gx, gx, z, theta = theta, phi = phi, zlim = zl, xlab = xlab,
                          ylab = ylab, zlab = "Valence", main = main, col = "grey95",
                          border = "grey60", ticktype = "detailed", ...)
    ok <- !is.na(v)
    pt <- grDevices::trans3d(xy[ok, 1], xy[ok, 2], v[ok], pm)
    above <- v[ok] >= x$height[k, ok]
    graphics::points(pt, pch = ifelse(above, 16, 1))
    if (labels) {
      graphics::text(pt, labels = x$labels[ok], pos = 3, cex = 0.8)
    }
  }
  return(invisible(list(x = gx, y = gx, z = z)))
}
