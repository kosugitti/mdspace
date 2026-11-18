# フォンミーゼスモデル(AMISESCAL)の試作関数 --------------------------------------
# 20261118_vonmises_try.R と 20261118_vonmises_stan.R から source する
# 設計メモ: develop/設計メモ_フォンミーゼスモデル_20261118.md
# 関数は base R で書き，固まったら R/vonmises_mds.R へ移す
#
# モデル(2次元):
#   theta_ij = x_i から x_j を見た角度，delta_ij = ||x_i - x_j||
#   pi_ij    = c_i * exp{kappa_i (cos(theta_ij - mu_i) - 1)}
#              (荘島の omega_i exp{kappa cos()} / I0(kappa) を山の高さ c_i で書き直したもの。
#               c_i = omega_i e^kappa / I0(kappa)。c_i < 1 で pi < 1 が守られる)
#   xi_ij    = (1 - pi_ij) delta_ij
# ストレス: sum (g_ij - xi_ij)^2 - lambda * U
#   U = sum_i sum_j max(g_ji - g_ij, 0) * vM(theta_ij | mu_i, kappa_i)(荘島 2012。lambda = 0 で 2011 版)
# 母数: x_i (2), mu_i, a_i = logit(kappa_i / KMAX), b_i = logit c_i
# kappa に上限 KMAX を置く: U を入れると kappa -> Inf で目的関数がいくらでも下がるため
KMAX <- 20


# 母数のベクトルを分ける ----------------------------------------------------------
.vm_unpack <- function(par, n) {
  list(X = matrix(par[1:(2 * n)], n, 2),
       mu = par[2 * n + 1:n],
       kappa = KMAX * stats::plogis(par[3 * n + 1:n]),
       c = stats::plogis(par[4 * n + 1:n]))
}

# セルごとの量(内部) -------------------------------------------------------------
.vm_parts <- function(par, cells) {
  p <- .vm_unpack(par, cells$n)
  i <- cells$i
  j <- cells$j
  dx <- p$X[j, 1] - p$X[i, 1]
  dy <- p$X[j, 2] - p$X[i, 2]
  delta <- sqrt(dx^2 + dy^2)
  theta <- atan2(dy, dx)
  cs <- cos(theta - p$mu[i])
  sn <- sin(theta - p$mu[i])
  pi_ij <- p$c[i] * exp(p$kappa[i] * (cs - 1))
  xi <- (1 - pi_ij) * delta
  # U の中の通常の vM 密度(omega = 1/(2pi))。besselI(expon.scaled) で桁あふれを避ける
  vm <- exp(p$kappa[i] * (cs - 1)) / (2 * pi * besselI(p$kappa[i], 0, expon.scaled = TRUE))
  list(p = p, dx = dx, dy = dy, delta = delta, theta = theta, cs = cs, sn = sn,
       pi = pi_ij, xi = xi, vm = vm)
}

# 目的関数 -------------------------------------------------------------------------
.vm_obj <- function(par, cells, lambda) {
  q <- .vm_parts(par, cells)
  sum((cells$g - q$xi)^2) - lambda * sum(cells$love * q$vm)
}

# 解析的な勾配 ---------------------------------------------------------------------
.vm_grad <- function(par, cells, lambda) {
  n <- cells$n
  q <- .vm_parts(par, cells)
  p <- q$p
  i <- cells$i
  j <- cells$j
  k <- p$kappa[i]
  e2 <- 2 * (q$xi - cells$g)               # dS/dxi
  # xi = (1 - pi) delta の各母数での微分
  dxi_dpi <- -q$delta
  dpi_dth <- -q$pi * k * q$sn
  dpi_dmu <- q$pi * k * q$sn
  dk_da <- k * (1 - k / KMAX)               # d kappa / d a
  dpi_da <- q$pi * (q$cs - 1) * dk_da
  dpi_db <- q$pi * (1 - p$c[i])            # d/d logit c
  # U の項: vm = exp{k(cs-1)} / (2 pi I0e(k))
  w <- lambda * cells$love
  i1 <- besselI(k, 1, expon.scaled = TRUE) / besselI(k, 0, expon.scaled = TRUE)
  dvm_dth <- -q$vm * k * q$sn
  dvm_dmu <- q$vm * k * q$sn
  dvm_da <- q$vm * ((q$cs - 1) - (i1 - 1)) * dk_da   # d log I0e / dk = I1/I0 - 1
  # 角度と距離の x_j での微分(x_i ではその符号を反転)
  d2 <- pmax(q$delta^2, 1e-12)
  ddel_x <- ifelse(q$delta > 1e-12, q$dx / q$delta, 0)
  ddel_y <- ifelse(q$delta > 1e-12, q$dy / q$delta, 0)
  dth_x <- -q$dy / d2
  dth_y <- q$dx / d2
  gth <- e2 * dxi_dpi * dpi_dth - w * dvm_dth   # dObj/dtheta
  gdel <- e2 * (1 - q$pi)                        # dObj/ddelta
  gx <- gdel * ddel_x + gth * dth_x
  gy <- gdel * ddel_y + gth * dth_y
  idx <- factor(c(j, i), levels = seq_len(n))
  GX <- cbind(as.vector(tapply(c(gx, -gx), idx, sum, default = 0)),
              as.vector(tapply(c(gy, -gy), idx, sum, default = 0)))
  fi <- factor(i, levels = seq_len(n))
  gmu <- as.vector(tapply(e2 * dxi_dpi * dpi_dmu - w * dvm_dmu, fi, sum, default = 0))
  ga <- as.vector(tapply(e2 * dxi_dpi * dpi_da - w * dvm_da, fi, sum, default = 0))
  gb <- as.vector(tapply(e2 * dxi_dpi * dpi_db, fi, sum, default = 0))
  c(GX, gmu, ga, gb)
}

# 初期値: 対称部の Torgerson 解，mu は片思いの相手の向きの平均 ------------------------------
.vm_start <- function(G, cells, random = FALSE) {
  n <- cells$n
  if (random) {
    X <- matrix(stats::rnorm(2 * n), n, 2)
    mu <- stats::runif(n, -pi, pi)
    a <- stats::qlogis(stats::runif(n, 0.5, 4) / KMAX)
    b <- stats::qlogis(stats::runif(n, 0.1, 0.6))
  } else {
    D <- G
    D[is.na(D)] <- t(D)[is.na(D)]
    D[is.na(D)] <- mean(D, na.rm = TRUE)
    diag(D) <- 0
    X <- stats::cmdscale((D + t(D)) / 2, k = 2)
    L <- pmax(t(D) - D, 0)                  # L_ij = i から j への片思いの量
    mu <- vapply(seq_len(n), function(r) {
      v <- colSums(L[r, ] * (t(X) - X[r, ]) / pmax(sqrt(colSums((t(X) - X[r, ])^2)), 1e-8))
      if (sum(abs(v)) < 1e-8) 0 else atan2(v[2], v[1])
    }, numeric(1))
    a <- rep(stats::qlogis(1 / KMAX), n)
    b <- rep(stats::qlogis(0.3), n)
  }
  c(X, mu, a, b)
}

# 識別のための正規化: 重心を原点，主軸を座標軸に回し，mu も同じだけ回す ----------------------
.vm_normalize <- function(p) {
  X <- sweep(p$X, 2, colMeans(p$X))
  V <- svd(X)$v
  if (det(V) < 0) V[, 2] <- -V[, 2]       # 回転だけ(鏡映はしない)
  rot <- atan2(V[2, 1], V[1, 1])
  p$X <- X %*% V
  p$mu <- ((p$mu - rot + pi) %% (2 * pi)) - pi
  p
}

# 本体 ----------------------------------------------------------------------------
vonmises_mds_try <- function(G, lambda = 0, nrandom = 20, seed = 1, maxit = 2000) {
  G <- as.matrix(G)
  n <- nrow(G)
  diag(G) <- NA
  # 尺度をそろえる: 非対角の二乗平均を 1 に(布置の単位がデータの単位になる)
  sc <- sqrt(mean(G^2, na.rm = TRUE))
  Gs <- G / sc
  use <- !is.na(Gs)
  cells <- list(n = n, i = row(Gs)[use], j = col(Gs)[use], g = Gs[use])
  cells$love <- pmax(t(Gs)[use] - Gs[use], 0)   # max(g_ji - g_ij, 0)
  cells$love[is.na(cells$love)] <- 0
  set.seed(seed)
  inits <- c(list(.vm_start(Gs, cells)),
             lapply(seq_len(nrandom), function(r) .vm_start(Gs, cells, random = TRUE)))
  fits <- lapply(inits, function(p0) {
    stats::optim(p0, .vm_obj, .vm_grad, cells = cells, lambda = lambda,
                 method = "BFGS", control = list(maxit = maxit, reltol = 1e-12))
  })
  vals <- vapply(fits, `[[`, numeric(1), "value")
  best <- fits[[which.min(vals)]]
  p <- .vm_normalize(.vm_unpack(best$par, n))
  q <- .vm_parts(best$par, cells)
  lab <- rownames(G)
  if (is.null(lab)) lab <- LETTERS[seq_len(n)]
  rownames(p$X) <- lab
  Pi <- Xi <- matrix(NA_real_, n, n, dimnames = list(lab, lab))
  Pi[use] <- q$pi
  Xi[use] <- q$xi * sc
  list(conf = p$X, mu = stats::setNames(p$mu, lab), kappa = stats::setNames(p$kappa, lab),
       peak = stats::setNames(p$c, lab), pi = Pi, xi = Xi, scale = sc,
       stress = sqrt(sum((cells$g - q$xi)^2) / sum(cells$g^2)),   # ストレス式1型
       objective = best$value, starts = vals, best_start = which.min(vals),
       converged = best$convergence == 0, lambda = lambda)
}

# 図: 点を中心に円を描き，その外側に pi_i(theta) を極座標で足す(荘島のスライドの形) -----------
plot_vm_try <- function(fit, main = "", rho = NULL) {
  X <- fit$conf
  if (is.null(rho)) rho <- 0.08 * diff(range(X))
  lim <- range(X) + c(-1, 1) * 3 * rho
  plot(X, asp = 1, xlim = lim, ylim = lim, type = "n", xlab = "Dim 1", ylab = "Dim 2", main = main)
  abline(h = 0, v = 0, col = "grey80", lty = 3)
  th <- seq(-pi, pi, length.out = 361)
  for (r in seq_len(nrow(X))) {
    lines(X[r, 1] + rho * cos(th), X[r, 2] + rho * sin(th), col = "grey50")
    dens <- fit$peak[r] * exp(fit$kappa[r] * (cos(th - fit$mu[r]) - 1))
    rr <- rho * (1 + 2 * dens)
    lines(X[r, 1] + rr * cos(th), X[r, 2] + rr * sin(th), col = "firebrick")
    segments(X[r, 1], X[r, 2], X[r, 1] + rho * (1 + 2 * fit$peak[r]) * cos(fit$mu[r]),
             X[r, 2] + rho * (1 + 2 * fit$peak[r]) * sin(fit$mu[r]), col = "firebrick")
  }
  text(X, labels = rownames(X), cex = 0.9)
}

