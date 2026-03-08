#' Dual Sacling with Mutual Averages
#'
#' @param dat A matrix or data frame
#' @param eps A criteria whether it reached conversion. Default value is 10e-5
#' @param iter.max A muximum numbers of iterations
#' @export


DualScalingMA <- function(dat, eps = 1e-10, iter.max = 1000) {
  nc <- ncol(dat)
  nr <- nrow(dat)
  cname <- colnames(dat)
  rname <- rownames(dat)
  mg.col <- colSums(dat)
  mg.row <- rowSums(dat)
  mg.t <- sum(mg.row)
  dm <- min(nc, nr) - 1
  normed.col <- matrix(ncol = nc, nrow = dm)
  projected.col <- matrix(ncol = nc, nrow = dm)
  normed.row <- matrix(ncol = nr, nrow = dm)
  projected.row <- matrix(ncol = nr, nrow = dm)
  singular <- c()

  dat.tmp <- dat
  # 0-exp
  for (i in 1:nrow(dat)) {
    for (j in 1:ncol(dat)) {
      dat.tmp[i, j] <- dat.tmp[i, j] - (mg.row[i] * mg.col[j] / mg.t)
    }
  }

  # 1-n
  for (i in 1:dm) {
    FLG <- FALSE
    iter <- 0
    u <- stats::rnorm(nc)
    tmp <- 0
    while (FLG == FALSE) {
      # algorithm
      v <- u %*% t(dat.tmp)
      v <- v / mg.row
      av <- as.vector((mg.row %*% t(v)) / mg.t)
      # v <- v - av * rep(1, nr)
      gy <- max(abs(v))
      v <- v / gy

      u <- v %*% dat.tmp
      u <- u / mg.col
      av <- as.vector((mg.col %*% t(u)) / mg.t)
      # u <- u - (av * rep(1, nc))
      gx <- max(abs(u))
      u <- u / gx

      # conv. check
      if (abs(tmp - gx) < eps) {
        FLG <- TRUE
      } else {
        tmp <- gx
      }
      if (iter > iter.max) {
        FLG <- TRUE
      } else {
        iter <- iter + 1
      }
    }


    if (u[which.max(abs(u))] < 0) u * -1
    if (v[which.max(abs(v))] < 0) v * -1

    eta2 <- gx * gy # Correlation ratio
    nu.c <- sqrt(sum(mg.col) / mg.col %*% t(u^2)) # multiplied constant
    nv.c <- sqrt(sum(mg.row) / mg.row %*% t(v^2))
    nu <- nu.c %*% u # normed u
    nv <- nv.c %*% v # normed v
    pu <- nu * sqrt(eta2)
    pv <- nv * sqrt(eta2)

    if (iter > iter.max) {
      cat("Warning::did not converge.")
    } else {
      singular[i] <- sqrt(eta2)
      normed.col[i, ] <- nu
      normed.row[i, ] <- nv
      projected.col[i, ] <- pu
      projected.row[i, ] <- pv
    }

    # residual Matrix
    dat.tmp <- dat.tmp - (mg.row %*% t(mg.col) / mg.t) * ((sqrt(eta2) * t(nv) %*% nu))
  }

  # Correct dims
  delta_k <- singular^2 / sum(singular^2) * 100
  cum_delta_k <- c()
  cum_delta_k[1] <- delta_k[1]
  for (i in 2:dm) {
    cum_delta_k[i] <- cum_delta_k[i - 1] + delta_k[i]
  }

  dm.cr <- which.max(cum_delta_k)
  if (dm.cr < dm) {
    singular <- singular[-(dm.cr + 1:dm)]
    delta_k <- delta_k[-(dm.cr + 1:dm)]
    cum_delta_k <- cum_delta_k[-(dm.cr + 1:dm)]
    normed.col <- normed.col[-(dm.cr + 1:dm), ]
    normed.row <- normed.row[-(dm.cr + 1:dm), ]
    projected.row <- projected.row[-(dm.cr + 1:dm), ]
    projected.col <- projected.col[-(dm.cr + 1:dm), ]
  }

  # var names
  colnames(normed.col) <- cname
  colnames(projected.col) <- cname
  colnames(normed.row) <- rname
  colnames(projected.row) <- rname
  dimName <- paste0("Dim", 1:dm.cr)
  rownames(normed.row) <- dimName
  rownames(normed.col) <- dimName
  rownames(projected.row) <- dimName
  rownames(projected.col) <- dimName

  # Style the output
  res <- t(cbind(singular, singular^2, delta_k, cum_delta_k))
  colnames(res) <- dimName
  rownames(res) <- c("singular value", "eigen value", "contribution", "cumulative contribution")

  normed.col <- t(normed.col)
  normed.row <- t(normed.row)
  projected.col <- t(projected.col)
  projected.row <- t(projected.row)



  return(list(
    result = res,
    ndims = dm.cr,
    singularValue = singular,
    eigenValue = singular^2,
    delta = delta_k,
    cumulativeDelta = cum_delta_k,
    NormedCol = normed.col,
    NormedRow = normed.row,
    ProjectedCol = projected.col,
    ProjectedRow = projected.row
  ))
}
