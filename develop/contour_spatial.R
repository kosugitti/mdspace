# 等高線マップモデルの空間統計: Moran の I と空間自己回帰(SAR)モデル -----------------------
# contour_mds() の結果(布置と人 x 対象の評定)を受けて，人ごとに
#   (段階0) Moran の I・帰無の期待値・並べ替え検定の p
#   (段階1) SAR モデル z = rho W z + e の最尤推定(mu・rho・sigma・rho = 0 との尤度比)
# を返す。W は Abelson と同じ 1 / (1 + d^2)(対角は 0)。rowstd = TRUE で行の和を 1 にそろえる。
# 設計: develop/設計メモ_等高線マップモデル_20261128.md の続き(先生 11/28: 両方の W を見る)

# 重み行列 ---------------------------------------------------------------------
contour_W <- function(conf, rowstd = FALSE) {
  W <- 1 / (1 + as.matrix(stats::dist(conf))^2)
  diag(W) <- 0
  if (rowstd) {
    W <- W / rowSums(W)
  }
  return(W)
}

# 順列を全部作る(n が小さいとき) ---------------------------------------------------
.all_perms <- function(n) {
  if (n == 1) {
    return(matrix(1L, 1, 1))
  }
  P <- .all_perms(n - 1)
  out <- lapply(seq_len(n), function(i) cbind(i, matrix(setdiff(seq_len(n), i)[P], nrow(P))))
  return(do.call(rbind, out))
}

# Moran の I(1 人分) -------------------------------------------------------------
.moran_I <- function(z, W) {
  return(length(z) / sum(W) * sum(z * (W %*% z)) / sum(z^2))
}

# 段階0: Moran の I ---------------------------------------------------------------
# nperm: n が max_exact 以下なら全部の並べ方を数える(正確な p)，それより大きければ nperm 回の乱数
contour_moran <- function(fit, rowstd = FALSE, max_exact = 8, nperm = 9999, seed = 1) {
  R <- fit$ratings
  res <- lapply(seq_len(nrow(R)), function(i) {
    ok <- !is.na(R[i, ])
    v <- R[i, ok]
    n <- length(v)
    W <- contour_W(fit$conf[ok, , drop = FALSE], rowstd)
    z <- v - mean(v)
    if (sum(z^2) == 0) {
      return(data.frame(person = rownames(R)[i], n = n, I = NA, E0 = -1 / (n - 1),
                        p_upper = NA, p_lower = NA, exact = NA))
    }
    I <- .moran_I(z, W)
    if (n <= max_exact) {
      P <- .all_perms(n)
      exact <- TRUE
    } else {
      set.seed(seed)
      P <- t(replicate(nperm, sample.int(n)))
      exact <- FALSE
    }
    Ip <- apply(P, 1, function(p) .moran_I(z[p], W))
    tol <- 1e-10
    data.frame(person = rownames(R)[i], n = n, I = I, E0 = -1 / (n - 1),
               p_upper = mean(Ip >= I - tol), p_lower = mean(Ip <= I + tol), exact = exact)
  })
  return(do.call(rbind, res))
}

# SAR の対数尤度(rho を与えると mu と sigma^2 は閉じた形で最大化) -------------------------
.sar_profile <- function(rho, v, W) {
  n <- length(v)
  A <- diag(n) - rho * W
  Av <- A %*% v
  X <- A %*% rep(1, n)
  mu <- sum(X * Av) / sum(X * X)
  e <- Av - X * mu
  s2 <- sum(e^2) / n
  ll <- as.numeric(determinant(A)$modulus) - n / 2 * (log(2 * pi * s2) + 1)
  return(list(ll = ll, mu = mu, sigma = sqrt(s2)))
}

# 段階1: SAR の最尤推定 --------------------------------------------------------------
contour_sar <- function(fit, rowstd = TRUE) {
  R <- fit$ratings
  res <- lapply(seq_len(nrow(R)), function(i) {
    ok <- !is.na(R[i, ])
    v <- R[i, ok]
    n <- length(v)
    W <- contour_W(fit$conf[ok, , drop = FALSE], rowstd)
    ev <- Re(eigen(W, only.values = TRUE)$values)
    lo <- 1 / min(ev) + 1e-8
    hi <- 1 / max(ev) - 1e-8
    if (stats::var(v) == 0) {
      return(data.frame(person = rownames(R)[i], n = n, mu = v[1], rho = NA, sigma = 0,
                        rho_lo = lo, rho_hi = hi, loglik = NA, LR = NA, p = NA))
    }
    o <- stats::optimize(function(r) -.sar_profile(r, v, W)$ll, c(lo, hi))
    best <- .sar_profile(o$minimum, v, W)
    null <- .sar_profile(0, v, W)
    LR <- 2 * (best$ll - null$ll)
    data.frame(person = rownames(R)[i], n = n, mu = best$mu, rho = o$minimum, sigma = best$sigma,
               rho_lo = lo, rho_hi = hi, loglik = best$ll, LR = LR,
               p = stats::pchisq(LR, 1, lower.tail = FALSE))
  })
  return(do.call(rbind, res))
}

# SAR から評定を作る(シミュレーション用) -----------------------------------------------
sar_simulate <- function(W, rho, mu = 0, sigma = 1) {
  n <- nrow(W)
  e <- stats::rnorm(n, 0, sigma)
  return(as.vector(mu + solve(diag(n) - rho * W, e)))
}
