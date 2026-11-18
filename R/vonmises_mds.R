# フォンミーゼスモデル(asymmetric von Mises scaling, AMISESCAL): 荘島 (2011, 2012) ----------

#' Asymmetric MDS by the von Mises model (AMISESCAL)
#'
#' Each object `i` is a point `x_i` in two dimensions with a von Mises
#' distribution of directions, with mean direction `mu_i` and concentration
#' `kappa_i`, which expresses toward which direction the object's one-sided
#' relations point (Shojima, 2011, 2012). With `theta_ij` the direction from
#' `x_i` to `x_j` and `delta_ij` their Euclidean distance, the model for the
#' dissimilarity from `i` to `j` is
#' `xi_ij = (1 - pi_ij) delta_ij`, where
#' `pi_ij = c_i exp{kappa_i (cos(theta_ij - mu_i) - 1)}`.
#' `pi_ij` is the von Mises density `omega_i exp{kappa_i cos(theta - mu_i)} / I0(kappa_i)`
#' of Shojima (2012) rewritten with its peak `c_i = omega_i exp(kappa_i) / I0(kappa_i)`;
#' `c_i < 1` keeps `pi_ij < 1`. An object `j` lying in the direction `mu_i`
#' looks closer from `i` than `i` looks from `j`.
#'
#' With `omega = "fixed"` the density is the ordinary von Mises density
#' (`omega_i = 1 / (2 pi)`, Shojima, 2011), so that `c_i` is determined by
#' `kappa_i`, and `kappa_i < 6.5` keeps `c_i < 1`. With `omega = "free"`
#' (Shojima, 2012) `c_i` is a parameter of its own. Leaving it free without
#' further constraint lets the solution degenerate (points move apart while
#' all `pi_ij` approach 1), which the original program prevents by an
#' additional term in the stress function. Here the estimates are instead the
#' mode of a posterior distribution (maximum a posteriori, MAP) with weak
#' priors, the same as the Bayesian version written in Stan
#' (`develop/vonmises_mds.stan` in the source repository):
#' `g_ij ~ N(xi_ij, sigma^2)`, deviations of the coordinates from their
#' centroid `~ N(0, prior$x_sd^2)`,
#' `kappa_i ~ Gamma(prior$kappa_shape, prior$kappa_rate)` truncated at the
#' upper bound, `c_i ~ Beta(prior$c_a, prior$c_b)` and
#' `sigma ~ half-t(4, 0, 1)`. The dissimilarities are divided by their root
#' mean square before fitting, so that the priors are on a common scale.
#'
#' The configuration is identified as in the Bayesian version by a Cholesky
#' factor (lower triangular with positive diagonal; Stan Development Team,
#' Reference Manual, "Cholesky factors of covariance matrices"): one anchor
#' object is placed at the origin, the second on the positive side of the
#' first axis and the third on the positive side of the second axis. The
#' anchors are chosen from the Torgerson solution of the symmetric part so that
#' they are far apart and not on a line. This anchoring is only a device for
#' computation: the likelihood and the prior (placed on the deviations from
#' the centroid, which is isotropic) do not change under translation,
#' rotation and reflection, so the solution does not depend on which objects
#' are the anchors. The mean directions rotate with the
#' configuration, so they are identified too. For presentation the returned
#' configuration is centered and rotated to its principal axes, and the mean
#' directions are rotated with it.
#'
#' The posterior mode has many local optima, so the rational start (the
#' Torgerson solution with each mean direction pointing at the objects toward
#' which the object has one-sided relations) and `nrandom` random starts are
#' tried and the best is returned. The minimization uses the analytic
#' gradient and [stats::optim()] (BFGS). With `hessian = TRUE` approximate
#' standard errors of `mu`, `kappa` and `c` are obtained from the Hessian of
#' the negative log posterior (Laplace approximation; delta method for the
#' transformed parameters).
#'
#' @param x A square matrix of asymmetric proximities. Rows evaluate columns
#'   (`x[i, j]` is from `i` to `j`). The diagonal is ignored; other missing
#'   cells are left out of the fit.
#' @param type `"dissimilarity"` (default) or `"similarity"`. Similarities
#'   are reversed as `max + min - x`. If `x` has a `type` attribute set by
#'   [as_proximity()], that is used.
#' @param omega `"free"` (default; Shojima, 2012) or `"fixed"` (`1 / (2 pi)`;
#'   Shojima, 2011). See Details.
#' @param prior A list with `x_sd`, `kappa_shape`, `kappa_rate`, `c_a` and
#'   `c_b`. See Details. Missing elements take the defaults
#'   `list(x_sd = 1, kappa_shape = 2, kappa_rate = 0.5, c_a = 2, c_b = 4)`.
#' @param nrandom Number of random starts in addition to the rational start.
#' @param maxit,reltol Passed to [stats::optim()].
#' @param hessian Logical; compute approximate standard errors.
#' @param seed Random seed for the random starts. The state of the random
#'   number generator is restored afterwards. `NULL` uses (and advances) the
#'   current state.
#' @return An object of class `c("mdspace_vonmises", "mdspace")` with `conf`
#'   (n x 2, in units of the scaled dissimilarities), `mu` (radians),
#'   `kappa`, `c` (peaks of `pi`), `sigma`, `stress`
#'   (`sqrt(sum((g - xi)^2) / sum(g^2))` on the scaled dissimilarities),
#'   `logpost` (log posterior at the mode, up to a constant), `niter`,
#'   `pi` and `xi` (n x n; `xi` in the original units), `scale` (the root mean
#'   square used for scaling), `anchor` (labels of the three anchors),
#'   `starts` (`logpost` of each start, the rational one first),
#'   `best_start`, `converged`, and with `hessian = TRUE`, `se` (a data frame
#'   of standard errors of `mu`, `kappa` and `c`).
#' @references Shojima, K. (2011). Asymmetric von Mises scaling.
#'   *Proceedings of the 39th Annual Meeting of the Behaviormetric Society of
#'   Japan*, 261-262.
#'
#' Shojima, K. (2012). On the stress function of asymmetric von Mises
#'   scaling. Paper presented at the Joint Meeting of the German
#'   Classification Society and the Japanese Classification Society.
#'
#' Mardia, K. V., & Jupp, P. E. (2000). *Directional statistics*. Wiley.
#'
#' Stan Development Team. *Stan reference manual*: Cholesky factors of
#'   covariance matrices. <https://mc-stan.org/docs/reference-manual/>
#' @examples
#' x <- matrix(c(NA, 5, 4, 5, 3,
#'               6, NA, 4, 2, 1,
#'               4, 4, NA, 3, 3,
#'               4, 1, 2, NA, 1,
#'               7, 1, 2, 1, NA), 5, byrow = TRUE,
#'             dimnames = list(LETTERS[1:5], LETTERS[1:5]))
#' fit <- vonmises_mds(x, nrandom = 2)
#' fit
#' @export
vonmises_mds <- function(x, type = c("dissimilarity", "similarity"),
                         omega = c("free", "fixed"), prior = list(),
                         nrandom = 10, maxit = 2000, reltol = 1e-10,
                         hessian = FALSE, seed = 1) {
  cl <- match.call()
  if (!is.null(attr(x, "type")) && missing(type)) {
    type <- attr(x, "type")
  }
  type <- match.arg(type)
  omega <- match.arg(omega)
  pr <- list(x_sd = 1, kappa_shape = 2, kappa_rate = 0.5, c_a = 2, c_b = 4)
  pr[names(prior)] <- prior
  X0 <- as_proximity(x, type = type)
  n <- nrow(X0)
  if (n < 3) {
    stop("At least 3 objects are needed.")
  }
  G <- unclass(X0)
  attr(G, "type") <- NULL
  diag(G) <- NA
  if (type == "similarity") {
    G <- max(G, na.rm = TRUE) + min(G, na.rm = TRUE) - G
  }
  if (any(G < 0, na.rm = TRUE)) {
    stop("Dissimilarities must be nonnegative.")
  }

  # 基準の3点を選び，先頭に並べ替える -------------------------------------------------
  ord <- .vm_anchor(G)
  Gs <- G[ord, ord]
  sc <- sqrt(mean(Gs^2, na.rm = TRUE))
  Gs <- Gs / sc
  use <- !is.na(Gs)
  if (sum(use) < 3) {
    stop("Too few observed off-diagonal cells.")
  }
  free <- omega == "free"
  cells <- list(n = n, i = row(Gs)[use], j = col(Gs)[use], g = Gs[use],
                free = free, kmax = if (free) 20 else 6.5, prior = pr)

  # 初期値: 合理的な初期値と乱数の初期値 ----------------------------------------------
  starts <- list(.vm_rational(Gs, cells))
  if (nrandom > 0) {
    if (!is.null(seed)) {
      # 呼び出し側の乱数の状態を壊さないよう，終わったら元に戻す
      old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
      on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()), add = TRUE)
      set.seed(seed)
    }
    for (r in seq_len(nrandom)) {
      starts[[r + 1]] <- .vm_random(cells)
    }
  }
  fits <- lapply(starts, function(p0) {
    stats::optim(p0, .vm_nlp, .vm_gr, cells = cells, method = "BFGS",
                 control = list(maxit = maxit, reltol = reltol))
  })
  vals <- vapply(fits, function(f) f$value, numeric(1))
  best_start <- which.min(vals)
  best <- fits[[best_start]]

  # 結果を組み立てる(元の並びに戻し，主軸に回す) -------------------------------------------
  p <- .vm_unpack(best$par, cells)
  q <- .vm_cells(p, cells)
  back <- order(ord)
  lab <- rownames(X0)
  Xc <- sweep(p$X, 2, colMeans(p$X))
  V <- svd(Xc)$v
  if (det(V) < 0) {
    V[, 2] <- -V[, 2]
  }
  rot <- atan2(V[2, 1], V[1, 1])
  conf <- (Xc %*% V)[back, , drop = FALSE]
  dimnames(conf) <- list(lab, c("D1", "D2"))
  mu <- (p$mu - rot + pi) %% (2 * pi) - pi
  Pi <- Xi <- matrix(NA_real_, n, n)
  Pi[use] <- q$pi
  Xi[use] <- q$xi * sc
  Pi <- Pi[back, back]
  Xi <- Xi[back, back]
  dimnames(Pi) <- dimnames(Xi) <- list(lab, lab)

  se <- NULL
  if (hessian) {
    se <- .vm_se(best$par, cells, p)[back, , drop = FALSE]
    rownames(se) <- lab
  }
  out <- .new_mdspace("vonmises", conf = conf,
                      stress = sqrt(sum((cells$g - q$xi)^2) / sum(cells$g^2)),
                      niter = unname(best$counts["gradient"]), data = X0, call = cl,
                      ndim = 2,
                      mu = stats::setNames(mu[back], lab),
                      kappa = stats::setNames(p$kappa[back], lab),
                      c = stats::setNames(p$c[back], lab),
                      sigma = p$sigma, logpost = -best$value,
                      pi = Pi, xi = Xi, scale = sc, omega = omega, prior = pr,
                      anchor = lab[ord[1:3]], starts = -vals, best_start = best_start,
                      converged = best$convergence == 0, se = se, type = type)
  return(out)
}

# 基準の3点(内部) ---------------------------------------------------------------------
# 対称部の Torgerson 解で，重心から最も遠い点→そこから最も遠い点→2点を結ぶ直線から最も遠い点
.vm_anchor <- function(G) {
  D <- .vm_fill(G)
  Z <- .vm_torgerson(D)
  a1 <- which.max(rowSums(sweep(Z, 2, colMeans(Z))^2))
  a2 <- which.max(rowSums(sweep(Z, 2, Z[a1, ])^2))
  u <- Z[a2, ] - Z[a1, ]
  W <- sweep(Z, 2, Z[a1, ])
  area <- abs(W[, 1] * u[2] - W[, 2] * u[1])
  area[c(a1, a2)] <- -Inf
  a3 <- which.max(area)
  return(c(a1, a2, a3, setdiff(seq_len(nrow(G)), c(a1, a2, a3))))
}

# 欠測を転置の値，それもなければ平均で埋める(初期値と基準点の選択だけに使う)(内部) ------------------
.vm_fill <- function(G) {
  D <- G
  D[is.na(D)] <- t(D)[is.na(D)]
  D[is.na(D)] <- mean(D, na.rm = TRUE)
  diag(D) <- 0
  return(D)
}

.vm_torgerson <- function(D) {
  Z <- suppressWarnings(stats::cmdscale((D + t(D)) / 2, k = 2))
  Z <- matrix(Z, nrow = nrow(D))
  if (ncol(Z) < 2) {
    Z <- cbind(Z, 0)
  }
  return(Z)
}

# 母数のベクトル(内部) -------------------------------------------------------------------
# u (2n - 3): 点2 = (exp(u1), 0)，点3 = (u2, exp(u3))，点4以降は自由。点1は原点
# mu (n)，a (n): kappa = kmax * plogis(a)，b (n, omega 自由のときだけ): c = plogis(b)，s: sigma = exp(s)
.vm_unpack <- function(par, cells) {
  n <- cells$n
  k <- 2 * n - 3
  u <- par[seq_len(k)]
  X <- matrix(0, n, 2)
  X[2, 1] <- exp(u[1])
  X[3, 1] <- u[2]
  X[3, 2] <- exp(u[3])
  if (n > 3) {
    X[4:n, ] <- matrix(u[-(1:3)], n - 3, 2, byrow = TRUE)
  }
  mu <- par[k + seq_len(n)]
  kappa <- cells$kmax * stats::plogis(par[k + n + seq_len(n)])
  if (cells$free) {
    c <- stats::plogis(par[k + 2 * n + seq_len(n)])
  } else {
    c <- 1 / (2 * pi * besselI(kappa, 0, expon.scaled = TRUE))
  }
  sigma <- exp(par[length(par)])
  return(list(X = X, mu = mu, kappa = kappa, c = c, sigma = sigma))
}

# 布置などから母数のベクトルを作る(内部) -------------------------------------------------
.vm_pack <- function(X, mu, kappa, c, sigma, cells) {
  n <- cells$n
  # コレスキーの枠へ: 点1を原点，点2を +x 軸，点3を上半面に
  X <- sweep(X, 2, X[1, ])
  ang <- atan2(X[2, 2], X[2, 1])
  # 角度 -ang の回転(行ベクトルに右から掛ける)
  R <- matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), 2)
  X <- X %*% R
  mu <- mu - ang
  if (X[3, 2] < 0) {
    X[, 2] <- -X[, 2]
    mu <- -mu
  }
  u <- c(log(max(X[2, 1], 1e-6)), X[3, 1], log(max(X[3, 2], 1e-6)))
  if (n > 3) {
    u <- c(u, as.vector(t(X[4:n, , drop = FALSE])))
  }
  kappa <- pmin(pmax(kappa, 1e-3), cells$kmax * 0.99)
  par <- c(u, mu, stats::qlogis(kappa / cells$kmax))
  if (cells$free) {
    par <- c(par, stats::qlogis(pmin(pmax(c, 0.01), 0.99)))
  }
  return(c(par, log(sigma)))
}

# 合理的な初期値: Torgerson 解，mu は片思いの相手の向きの重みつき平均(内部) --------------------
.vm_rational <- function(Gs, cells) {
  D <- .vm_fill(Gs)
  X <- .vm_torgerson(D)
  L <- pmax(t(D) - D, 0)
  mu <- vapply(seq_len(cells$n), function(r) {
    W <- t(X) - X[r, ]
    len <- pmax(sqrt(colSums(W^2)), 1e-8)
    v <- colSums(L[r, ] * t(W) / len)
    if (sum(abs(v)) < 1e-8) 0 else atan2(v[2], v[1])
  }, numeric(1))
  return(.vm_pack(X, mu, rep(1, cells$n), rep(0.3, cells$n), 0.5, cells))
}

.vm_random <- function(cells) {
  n <- cells$n
  X <- matrix(stats::rnorm(2 * n), n, 2)
  return(.vm_pack(X, stats::runif(n, -pi, pi), stats::runif(n, 0.5, 4),
                  stats::runif(n, 0.1, 0.6), 0.5, cells))
}

# セルごとの量(内部) -------------------------------------------------------------------
.vm_cells <- function(p, cells) {
  i <- cells$i
  j <- cells$j
  dx <- p$X[j, 1] - p$X[i, 1]
  dy <- p$X[j, 2] - p$X[i, 2]
  delta <- sqrt(dx^2 + dy^2)
  theta <- atan2(dy, dx)
  cs <- cos(theta - p$mu[i])
  sn <- sin(theta - p$mu[i])
  pi_ij <- p$c[i] * exp(p$kappa[i] * (cs - 1))
  return(list(dx = dx, dy = dy, delta = delta, cs = cs, sn = sn, pi = pi_ij,
              xi = (1 - pi_ij) * delta))
}

# 負の対数事後密度(定数を除く・ヤコビアンなし＝Stan の optimize と同じ)(内部) ----------------------
.vm_nlp <- function(par, cells) {
  p <- .vm_unpack(par, cells)
  q <- .vm_cells(p, cells)
  pr <- cells$prior
  S <- sum((cells$g - q$xi)^2)
  L <- length(cells$g)
  out <- L * log(p$sigma) + S / (2 * p$sigma^2) +
    sum(sweep(p$X, 2, colMeans(p$X))^2) / (2 * pr$x_sd^2) -
    sum((pr$kappa_shape - 1) * log(p$kappa) - pr$kappa_rate * p$kappa) +
    2.5 * log(1 + p$sigma^2 / 4)
  if (cells$free) {
    out <- out - sum((pr$c_a - 1) * log(p$c) + (pr$c_b - 1) * log(1 - p$c))
  }
  return(out)
}

# その勾配(内部) -----------------------------------------------------------------------
.vm_gr <- function(par, cells) {
  n <- cells$n
  p <- .vm_unpack(par, cells)
  q <- .vm_cells(p, cells)
  pr <- cells$prior
  i <- cells$i
  j <- cells$j
  sig2 <- p$sigma^2
  e <- q$xi - cells$g
  w <- e / sig2                                  # d NLP / d xi
  k <- p$kappa[i]
  # 布置: xi は x_j - x_i だけを通して布置に依存する
  d2 <- pmax(q$delta^2, 1e-12)
  ok <- q$delta > 1e-12
  gdel <- w * (1 - q$pi)
  gth <- w * (-q$delta) * (-q$pi * k * q$sn)
  gx <- ifelse(ok, gdel * q$dx / q$delta, 0) + gth * (-q$dy / d2)
  gy <- ifelse(ok, gdel * q$dy / q$delta, 0) + gth * (q$dx / d2)
  idx <- factor(c(j, i), levels = seq_len(n))
  GX <- cbind(as.vector(tapply(c(gx, -gx), idx, sum, default = 0)),
              as.vector(tapply(c(gy, -gy), idx, sum, default = 0)))
  # 事前分布は重心からのずれに置く(平行移動で変わらない)。d/dX sum((X - mean)^2) / 2 = X - mean
  GX <- GX + sweep(p$X, 2, colMeans(p$X)) / pr$x_sd^2
  gu <- c(GX[2, 1] * p$X[2, 1], GX[3, 1], GX[3, 2] * p$X[3, 2])
  if (n > 3) {
    gu <- c(gu, as.vector(t(GX[4:n, , drop = FALSE])))
  }
  fi <- factor(i, levels = seq_len(n))
  wp <- w * (-q$delta)                           # d NLP / d pi
  gmu <- as.vector(tapply(wp * q$pi * k * q$sn, fi, sum, default = 0))
  # kappa: pi = c exp{kappa (cs - 1)}。omega 固定なら c も kappa の関数(dc/dkappa = c (1 - I1/I0))
  dpi_dk <- q$pi * (q$cs - 1)
  if (!cells$free) {
    A <- besselI(k, 1, expon.scaled = TRUE) / besselI(k, 0, expon.scaled = TRUE)
    dpi_dk <- dpi_dk + q$pi * (1 - A)
  }
  gk <- as.vector(tapply(wp * dpi_dk, fi, sum, default = 0)) -
    ((pr$kappa_shape - 1) / p$kappa - pr$kappa_rate)
  ga <- gk * p$kappa * (1 - p$kappa / cells$kmax)
  out <- c(gu, gmu, ga)
  if (cells$free) {
    gc <- as.vector(tapply(wp * q$pi / p$c[i], fi, sum, default = 0)) -
      ((pr$c_a - 1) / p$c - (pr$c_b - 1) / (1 - p$c))
    out <- c(out, gc * p$c * (1 - p$c))
  }
  S <- sum(e^2)
  L <- length(e)
  gs <- L - S / sig2 + 2.5 * (p$sigma^2 / 2) / (1 + p$sigma^2 / 4)
  return(c(out, gs))
}

# ヘッセ行列からの近似標準誤差(内部) --------------------------------------------------------
.vm_se <- function(par, cells, p) {
  n <- cells$n
  H <- stats::optimHess(par, .vm_nlp, .vm_gr, cells = cells)
  V <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(V)) {
    warning("The Hessian is singular; standard errors are not available.")
    return(data.frame(mu = rep(NA_real_, n), kappa = NA_real_, c = NA_real_))
  }
  k <- 2 * n - 3
  v <- pmax(diag(V), 0)
  se_mu <- sqrt(v[k + seq_len(n)])
  dk <- p$kappa * (1 - p$kappa / cells$kmax)
  se_a <- sqrt(v[k + n + seq_len(n)])
  se_k <- dk * se_a
  if (cells$free) {
    se_c <- p$c * (1 - p$c) * sqrt(v[k + 2 * n + seq_len(n)])
  } else {
    A <- besselI(p$kappa, 1, expon.scaled = TRUE) / besselI(p$kappa, 0, expon.scaled = TRUE)
    se_c <- abs(p$c * (1 - A)) * se_k
  }
  return(data.frame(mu = se_mu, kappa = se_k, c = se_c))
}

#' @export
print.mdspace_vonmises <- function(x, digits = 4, ...) {
  NextMethod()
  cat("omega:", x$omega, if (x$omega == "fixed") "(1 / (2 pi); Shojima, 2011)" else "(free; Shojima, 2012)", "\n")
  cat("sigma:", format(round(x$sigma, digits), nsmall = digits),
      " log posterior:", format(round(x$logpost, 2), nsmall = 2), "\n")
  if (!isTRUE(x$converged)) {
    cat("Note: the optimizer did not report convergence.\n")
  }
  if (length(x$starts) > 1) {
    cat(sprintf("Starts: %d (rational + %d random); best = %s\n",
                length(x$starts), length(x$starts) - 1,
                if (x$best_start == 1) "rational" else paste("random", x$best_start - 1)))
  }
  return(invisible(x))
}

#' @export
summary.mdspace_vonmises <- function(object, ...) {
  out <- NextMethod()
  tab <- cbind(object$conf, mu_deg = object$mu * 180 / pi, kappa = object$kappa, c = object$c)
  if (!is.null(object$se)) {
    tab <- cbind(tab, se_mu_deg = object$se$mu * 180 / pi, se_kappa = object$se$kappa,
                 se_c = object$se$c)
  }
  out$conf <- tab
  return(out)
}

#' Plot a von Mises model solution
#'
#' Each object is drawn as a point with a small circle, and the curve outside
#' the circle shows `pi_i(theta)`, the strength of its one-sided relations in
#' each direction (as in the figures of Shojima, 2011, 2012). A line from the
#' point shows the mean direction `mu_i`.
#'
#' @param x An object of class `mdspace_vonmises`.
#' @param labels Logical; draw object labels.
#' @param rho Radius of the base circles. `NULL` uses 6% of the range of the
#'   configuration.
#' @param ... Passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
plot.mdspace_vonmises <- function(x, labels = TRUE, rho = NULL, ...) {
  X <- x$conf
  if (is.null(rho)) {
    rho <- 0.06 * diff(range(X))
  }
  lim <- range(X) + c(-1, 1) * 3.2 * rho
  graphics::plot(X, asp = 1, xlim = lim, ylim = lim, type = "n",
                 xlab = "Dimension 1", ylab = "Dimension 2",
                 main = "von Mises model", ...)
  graphics::abline(h = 0, v = 0, col = "grey80", lty = 3)
  th <- seq(-pi, pi, length.out = 361)
  for (r in seq_len(nrow(X))) {
    graphics::lines(X[r, 1] + rho * cos(th), X[r, 2] + rho * sin(th), col = "grey50")
    dens <- x$c[r] * exp(x$kappa[r] * (cos(th - x$mu[r]) - 1))
    rr <- rho * (1 + 2 * dens)
    graphics::lines(X[r, 1] + rr * cos(th), X[r, 2] + rr * sin(th), col = "firebrick")
    tip <- rho * (1 + 2 * x$c[r])
    graphics::segments(X[r, 1], X[r, 2], X[r, 1] + tip * cos(x$mu[r]),
                       X[r, 2] + tip * sin(x$mu[r]), col = "firebrick")
  }
  if (labels) {
    graphics::text(X, labels = x$labels, cex = 0.8)
  }
  return(invisible(x))
}
