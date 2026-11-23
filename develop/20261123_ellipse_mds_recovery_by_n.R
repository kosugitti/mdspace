# 楕円モデルの半軸は，対象がいくつあれば落ち着いて推定できるか --------------------------------------------
# 対象 12 個では楕円の形が stress でほとんど決まらない(細長く大きな楕円へ流れる)。対象 j の楕円の形は
# 「j から見た k の方向ごとの半径」で決まるので，対象が増えれば方向が増えて決まるはず，という見立てを確かめる。
# 楕円モデル(Okada 1990 の式 1〜3)から人工データを作り，対象の数 N と誤差 E を変えて ellipse_mds() の回復を見る
#   布置: 2 次元。次元の分散を変えて(SD 1.3 と 0.7)次元の順が入れ替わりにくくし，重心 0・平方和 N に正規化
#   半軸: log r_jt = c_j + e_jt，c_j ~ N(log 0.4, 0.4)(大きさ)，e_jt ~ N(0, 0.4) を対象ごとに平均 0 に(形)
#   誤差: SD = E * SD(m)(Okada & Imaizumi 1987 の E を SD の比として使う。20261122 の注と同じ)
# 回復の指標: 点間距離・大きさ(半軸の幾何平均の対数)・形(log r_j1 - log r_j2)の相関，半軸の最大値

pacman::p_load(tidyverse)
devtools::load_all(quiet = TRUE)

# 楕円モデルから非類似度を作る -----------------------------------------------------------------------
simulate_ellipse <- function(N, E, seed, P = 2) {
  set.seed(seed)
  Z <- matrix(rnorm(N * P), N, P) %*% diag(c(1.3, 0.7))
  Z <- sweep(Z, 2, colMeans(Z))
  Z <- Z * sqrt(N / sum(Z^2))
  size <- rnorm(N, log(0.4), 0.4)
  shape <- matrix(rnorm(N * P, 0, 0.4), N, P)
  shape <- shape - rowMeans(shape)
  r <- exp(size + shape)
  M <- matrix(0, N, N)
  for (j in seq_len(N)) {
    for (k in seq_len(N)) {
      if (j == k) next
      dif <- Z[j, ] - Z[k, ]
      d <- sqrt(sum(dif^2))
      M[j, k] <- d - d / sqrt(sum((dif / r[j, ])^2)) + d / sqrt(sum((dif / r[k, ])^2))
    }
  }
  off <- row(M) != col(M)
  d <- M + matrix(rnorm(N * N, 0, E * sd(M[off])), N, N)
  diag(d) <- NA
  dimnames(d) <- list(paste0("o", 1:N), paste0("o", 1:N))
  list(d = d, Z = Z, r = r, size = size, shape = shape[, 1] - shape[, 2])
}

# 当てはめて回復を測る(次元の入れ替えは座標の相関で合わせる) ----------------------------------------
fit_recover <- function(sim, nrandom = 3) {
  t0 <- proc.time()[["elapsed"]]
  out <- suppressWarnings(ellipse_mds(sim$d, ndim = 2, model = "II", type = "dissimilarity",
                                      nrandom = nrandom, check_orientation = FALSE))
  sec <- proc.time()[["elapsed"]] - t0
  cc <- abs(cor(out$conf, sim$Z))
  perm <- if (cc[1, 1] + cc[2, 2] >= cc[1, 2] + cc[2, 1]) 1:2 else 2:1   # 推定の次元 t が真の次元 perm[t]
  lr <- log(out$semiaxis)
  lr_true <- log(sim$r[, perm])
  tibble(stress = out$stress, sec = sec, swapped = perm[1] == 2,
         cor_dist = cor(c(dist(out$conf)), c(dist(sim$Z))),
         cor_size = cor(rowMeans(lr), rowMeans(lr_true)),
         cor_shape = cor(lr[, 1] - lr[, 2], lr_true[, 1] - lr_true[, 2]),
         cor_axis = cor(c(lr), c(lr_true)),
         aspect_max = max(exp(abs(lr[, 1] - lr[, 2]))),
         semiaxis_max = max(out$semiaxis),
         semiaxis_true_max = max(sim$r))
}

# 条件を振る(コマンドライン引数に time があれば，時間を測るだけ) ------------------------------------
if ("time" %in% commandArgs(trailingOnly = TRUE)) {
  res <- tibble(N = c(12, 30, 100)) %>%
    mutate(out = map(N, function(N) fit_recover(simulate_ellipse(N, 0.1, 1), nrandom = 1))) %>%
    unnest(out)
  saveRDS(res, "develop/20261123_ellipse_mds_recovery_by_n_time.rds")
} else {
  grid <- expand_grid(N = c(12, 20, 30, 50, 100), E = c(0.1, 0.3, 0.5), rep = 1:5) %>%
    mutate(seed = 3000 + row_number())
  res <- grid %>%
    mutate(out = pmap(list(N, E, seed), function(N, E, seed) fit_recover(simulate_ellipse(N, E, seed)))) %>%
    unnest(out)
  saveRDS(res, "develop/20261123_ellipse_mds_recovery_by_n.rds")
}
