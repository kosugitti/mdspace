# 二相三元の円モデル(Okada & Imaizumi, 1997)の内部関数 --------------------------------
# circle_mds() に行列のリストを渡したときに使う。
# m_jki = w_i d_jk - v_jki (r_j - r_k)，v_jki = w_i d_jk / sqrt(sum_t ((x_jt - x_kt) / u_it)^2)
# (nASYMS の出力で確認済み。w_i は m 全体に掛かる)
# S = sqrt(mean_i S_i^2)，S_i はソースごとのストレス式2。
# w_i は m 全体の定数倍なので S に効かない＝推定では w_i = 1 とし，最後に原典(3.7)の正規化で決める。
# 推定する母数は X (n x p)，r (n)，U (N x p)。u_it = exp(b_it - mean(b))，b_it = c_i + e_it とし，
# 全体の幾何平均を 1 に固定(U と r の大きさの不定性を止める)，ソースの中の形 e_it だけを
# [-L, L] に収める(縦横比 <= kappa＝楕円が線分につぶれる退化を止める)。

# ソース i の m と微分に要る量(w_i = 1) ----------------------------------------------------
.c3_msrc <- function(X, r, u, cl) {
  dif <- X[cl$j, , drop = FALSE] - X[cl$k, , drop = FALSE]
  d <- sqrt(rowSums(dif^2))
  q <- sqrt(rowSums(sweep(dif, 2, u, "/")^2))
  v <- ifelse(q > 1e-12, d / q, 0)
  return(list(m = d - v * (r[cl$j] - r[cl$k]), dif = dif, d = d, q = q, v = v))
}

.c3_unpack <- function(par, n, p, N) {
  return(list(X = matrix(par[seq_len(n * p)], n, p), r = par[n * p + seq_len(n)],
              U = matrix(par[n * p + n + seq_len(N * p)], N, p)))
}

# S と各ソースの S_i(内部) -------------------------------------------------------------
.c3_stress_parts <- function(par, n, p, N, cells) {
  P <- .c3_unpack(par, n, p, N)
  Si2 <- vapply(seq_len(N), function(i) {
    o <- .c3_msrc(P$X, P$r, P$U[i, ], cells[[i]])
    mh <- .monotone(o$m, cells[[i]]$s, cells[[i]]$decreasing)
    sum((o$m - mh)^2) / sum((o$m - mean(o$m))^2)
  }, numeric(1))
  return(Si2)
}

.c3_stress <- function(par, n, p, N, cells) {
  return(sqrt(mean(.c3_stress_parts(par, n, p, N, cells))))
}

# 勾配(X，r，U について) -------------------------------------------------------------------
.c3_grad <- function(par, n, p, N, cells) {
  P <- .c3_unpack(par, n, p, N)
  GX <- matrix(0, n, p)
  Gr <- numeric(n)
  GU <- matrix(0, N, p)
  Si2 <- numeric(N)
  parts <- vector("list", N)
  for (i in seq_len(N)) {
    cl <- cells[[i]]
    o <- .c3_msrc(P$X, P$r, P$U[i, ], cl)
    m <- o$m
    mh <- .monotone(m, cl$s, cl$decreasing)
    A <- sum((m - mh)^2)
    B <- sum((m - mean(m))^2)
    Si2[i] <- A / B
    parts[[i]] <- list(o = o, g = 2 * ((m - mh) - (A / B) * (m - mean(m))) / B)
  }
  S <- sqrt(mean(Si2))
  if (S < 1e-12) {
    return(numeric(length(par)))
  }
  for (i in seq_len(N)) {
    cl <- cells[[i]]
    o <- parts[[i]]$o
    g <- parts[[i]]$g / (2 * N * S)
    u <- P$U[i, ]
    dr <- P$r[cl$j] - P$r[cl$k]
    idx <- factor(c(cl$j, cl$k), levels = seq_len(n))
    ok <- o$d > 1e-12 & o$q > 1e-12
    for (t in seq_len(p)) {
      a <- o$dif[, t]
      dd <- ifelse(ok, a / o$d, 0)
      dq <- ifelse(ok, a / (u[t]^2 * o$q), 0)
      dv <- ifelse(ok, (dd * o$q - o$d * dq) / o$q^2, 0)
      dm <- dd - dv * dr
      GX[, t] <- GX[, t] + as.vector(tapply(c(g * dm, -g * dm), idx, sum, default = 0))
      dvu <- ifelse(ok, o$d * (a^2 / (u[t]^3 * o$q)) / o$q^2, 0)
      GU[i, t] <- sum(g * (-dvu * dr))
    }
    Gr <- Gr + as.vector(tapply(c(-g * o$v, g * o$v), idx, sum, default = 0))
  }
  return(c(GX, Gr, GU))
}

# 1つの初期値からの最小化(縦横比の制約つき・L-BFGS-B) --------------------------------------
.c3_optim <- function(X, r, cc, e, cells, n, N, kappa, maxit) {
  p <- ncol(X)
  L <- log(kappa) / 2
  ic <- n * p + n + seq_len(N)
  ie <- n * p + n + N + seq_len(N * p)
  iu <- n * p + n + seq_len(N * p)
  to_par <- function(th) {
    b <- th[ic] + th[ie]
    return(c(th[seq_len(n * p + n)], exp(b - mean(b))))
  }
  fn <- function(th) .c3_stress(to_par(th), n, p, N, cells)
  gr <- function(th) {
    par <- to_par(th)
    G <- .c3_grad(par, n, p, N, cells)
    gu <- G[iu] * par[iu]
    gb <- gu - mean(gu)
    return(c(G[seq_len(n * p + n)], rowSums(matrix(gb, N, p)), gb))
  }
  th0 <- c(X, r, cc, pmin(pmax(e, -L), L))
  lower <- c(rep(-Inf, n * p + n + N), rep(-L, N * p))
  upper <- c(rep(Inf, n * p + n + N), rep(L, N * p))
  res <- stats::optim(th0, fn, gr, method = "L-BFGS-B", lower = lower, upper = upper,
                      control = list(maxit = maxit, factr = 1e3))
  P <- .c3_unpack(to_par(res$par), n, p, N)
  e_hat <- matrix(res$par[ie], N, p)
  return(list(X = P$X, r = P$r, U = P$U, cc = res$par[ic], e = e_hat,
              stress = res$value, niter = unname(res$counts["gradient"]),
              converged = res$convergence == 0,
              at_bound = any(abs(abs(e_hat) - L) < 1e-6)))
}

# 初期値: ソースごとの非類似度(最小を 0 に)の平均を対称化し，加算定数つき Torgerson -------------
.c3_init <- function(Xs, p, type) {
  N <- length(Xs)
  deltas <- lapply(Xs, function(X) {
    D <- if (type == "similarity") -X else X
    diag(D) <- NA
    D <- D - min(D, na.rm = TRUE)
    miss <- is.na(D)
    D[miss] <- t(D)[miss]
    D[is.na(D)] <- mean(D, na.rm = TRUE)
    diag(D) <- 0
    D
  })
  delta <- Reduce(`+`, deltas) / N
  sym <- (delta + t(delta)) / 2
  skew <- (delta - t(delta)) / 2
  conf <- suppressWarnings(stats::cmdscale(sym, k = p, add = TRUE)$points)
  conf <- matrix(conf, nrow = nrow(delta))
  if (ncol(conf) < p) {
    conf <- cbind(conf, matrix(stats::rnorm(nrow(conf) * (p - ncol(conf)), 0, 0.01),
                               nrow(conf)))
  }
  radius <- colMeans(skew)
  return(list(conf = conf, radius = radius - min(radius)))
}

# 二相三元の本体(circle_mds から呼ぶ) -------------------------------------------------------
.circle3_fit <- function(Xs, ndim, type, ndim_max, nrandom, maxit, seed, kappa) {
  N <- length(Xs)
  n <- nrow(Xs[[1]])
  cells <- lapply(Xs, function(X) {
    use <- row(X) != col(X) & !is.na(X)
    list(s = X[use], j = row(X)[use], k = col(X)[use], decreasing = type == "similarity")
  })
  # 合理的初期値: 最大次元で解き，寄与(x の二乗和＋e の二乗和)の小さい次元を1つずつ落とす
  init <- .c3_init(Xs, ndim_max, type)
  X <- init$conf
  r <- init$radius
  cc <- rep(0, N)
  e <- matrix(0, N, ndim_max)
  path <- numeric(0)
  fit <- NULL
  for (p in seq(ndim_max, ndim)) {
    if (ncol(X) > p) {
      Xc <- sweep(X, 2, colMeans(X))
      contrib <- colSums(Xc^2) / sum(Xc^2) + colSums(e^2) / max(sum(e^2), 1e-12)
      keep <- sort(order(contrib, decreasing = TRUE)[seq_len(p)])
      X <- X[, keep, drop = FALSE]
      e <- e[, keep, drop = FALSE]
    }
    fit <- .c3_optim(X, r, cc, e, cells, n, N, kappa, maxit)
    X <- fit$X
    r <- fit$r
    cc <- fit$cc
    e <- fit$e
    path[as.character(p)] <- fit$stress
  }
  best <- fit
  best_start <- 1L
  starts <- fit$stress
  # 乱数の初期値 ----------------------------------------------------------------------
  if (nrandom > 0) {
    if (!is.null(seed)) {
      old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
      on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()), add = TRUE)
      set.seed(seed)
    }
    L <- log(kappa) / 2
    for (s in seq_len(nrandom)) {
      f <- .c3_optim(matrix(stats::rnorm(n * ndim), n, ndim), stats::runif(n, 0, 0.5),
                     stats::rnorm(N, 0, 0.3), matrix(stats::runif(N * ndim, -L, L), N, ndim),
                     cells, n, N, kappa, maxit)
      starts <- c(starts, f$stress)
      if (f$stress < best$stress) {
        best <- f
        best_start <- s + 1L
      }
    }
  }
  # 正規化(nASYMS と同じ形): 布置は重心 0・二乗和 n，半径は同じ倍率で最小 0 -----------------
  Xc <- sweep(best$X, 2, colMeans(best$X))
  a <- sqrt(n / sum(Xc^2))
  conf <- Xc * a
  radius <- (best$r - min(best$r)) * a
  # 対称重み w_i: m_jki = w_i (d_jk - (d_jk / q_jki)(r_j - r_k)) と w_i が全体に掛かるので，
  # ソースごとに sum (m - mean m)^2 = n となる倍率を w_i とする(原典 3.7・nASYMS の出力と同じ形)
  U <- best$U
  Si2 <- .c3_stress_parts(c(conf, radius, U), n, ndim, N, cells)
  w <- vapply(seq_len(N), function(i) {
    o <- .c3_msrc(conf, radius, U[i, ], cells[[i]])
    sqrt(n / sum((o$m - mean(o$m))^2))
  }, numeric(1))
  return(list(conf = conf, radius = radius, U = U, w = w, stress = best$stress,
              stress_source = sqrt(Si2), niter = best$niter, converged = best$converged,
              path = path, starts = starts, best_start = best_start,
              at_bound = best$at_bound, aspect = apply(U, 1, function(u) max(u) / min(u)),
              cells = cells))
}

# circle_mds(list) の入口(内部) ------------------------------------------------------------
.circle3 <- function(x, ndim, type, ndim_max, nrandom, maxit, seed, kappa, cl) {
  Xs <- lapply(x, function(m) {
    m <- as.matrix(m)
    storage.mode(m) <- "double"
    m
  })
  n <- nrow(Xs[[1]])
  if (any(vapply(Xs, function(m) nrow(m) != n || ncol(m) != n, logical(1)))) {
    stop("All matrices in `x` must be square and of the same size.")
  }
  if (length(Xs) < 2) {
    stop("A list of at least two matrices is needed for the two-mode three-way model.")
  }
  if (ndim < 1 || ndim >= n) {
    stop("`ndim` must be between 1 and n - 1.")
  }
  if (ndim_max < ndim || ndim_max >= n) {
    stop("`ndim_max` must be between `ndim` and n - 1.")
  }
  if (kappa < 1) {
    stop("`kappa` must be at least 1.")
  }
  lab <- rownames(Xs[[1]])
  if (is.null(lab)) lab <- colnames(Xs[[1]])
  if (is.null(lab)) lab <- as.character(seq_len(n))
  src <- names(x)
  if (is.null(src)) src <- paste0("S", seq_along(Xs))
  Xs <- lapply(Xs, function(m) {
    dimnames(m) <- list(lab, lab)
    m
  })
  fit <- .circle3_fit(Xs, ndim, type, ndim_max, nrandom, maxit, seed, kappa)
  if (fit$at_bound) {
    warning(sprintf("The aspect ratio of the ellipses reached the bound `kappa` = %g in some source(s); the data may favour collapsed ellipses.", kappa))
  }
  conf <- fit$conf
  dimnames(conf) <- list(lab, paste0("D", seq_len(ndim)))
  radius <- stats::setNames(fit$radius, lab)
  weights <- data.frame(source = src, w = fit$w, fit$U)
  names(weights) <- c("source", "w", paste0("u", seq_len(ndim)))
  out <- .new_mdspace("circle", conf = conf, stress = fit$stress, niter = fit$niter,
                      data = Xs, call = cl, ndim = ndim, radius = radius, type = type,
                      weights = weights, aspect = stats::setNames(fit$aspect, src),
                      stress_source = stats::setNames(fit$stress_source, src),
                      at_bound = fit$at_bound, path = fit$path, starts = fit$starts,
                      best_start = fit$best_start, converged = fit$converged,
                      kappa = kappa, threeway = TRUE)
  return(out)
}
