# 2. 千野の印字した全固有分解から各時点の行列を組み直し，gstatis_mds() の結果を照合する --------
#    (欠測の埋め方が特定できないので，千野が補完したあとの H を固有分解から復元して使う)
#    照合: 妥協の固有値・軌跡(千野は対象ごとに最大の長さで割って印字)・Lp ノルム・H と Xi の距離
dir <- "develop/references/chino1999_gstatis/"
ev <- read.csv(paste0(dir, "chino1999_eigenvalues.csv")); vec <- read.csv(paste0(dir, "chino1999_eigenvectors.csv"))
tr <- read.csv(paste0(dir, "chino1999_trajectory.csv")); lp <- read.csv(paste0(dir, "chino1999_lpnorm.csv")); dst <- read.csv(paste0(dir, "chino1999_distance.csv"))
rebuildX <- function(g, stage) { v <- vec[vec$group == g & vec$stage == stage & vec$type == "raw", ]; n <- max(v$object); H <- matrix(0i, n, n)
  for (r in unique(v$rank)) { u <- v[v$rank == r, ]; z <- complex(real = u$real, imaginary = u$imag); H <- H + u$eigenvalue[1] * z %*% Conj(t(z)) }
  X <- Re(H) + Im(H); dimnames(X) <- list(1:n, 1:n); X }
res <- NULL
for (g in 1:4) {
  times <- unique(ev$stage[ev$group == g & !(ev$stage %in% c("RV", "compromise"))])
  Xs <- setNames(lapply(times, function(tm) rebuildX(g, tm)), times)
  f <- gstatis_mds(Xs, ndim = 3, ncomp = 3)
  fneg <- gstatis_mds(lapply(Xs, function(m) -m), ndim = 3, ncomp = 3)  # 重みの符号を反転したのと同じ
  for (vp in 1:3) {
    # 固有ベクトルの符号を千野にそろえる(違えば H をすべて反転した解の同じ視点を使う)
    cw <- vec$real[vec$group == g & vec$stage == "RV" & vec$rank == vp]
    flip <- sum(cw * f$compromise[[vp]]$weights) < 0
    cp <- if (flip) fneg$compromise[[vp]] else f$compromise[[vp]]
    sg <- 1
    ce <- sort(ev$eigenvalue[ev$group == g & ev$stage == "compromise" & ev$viewpoint == vp], decreasing = TRUE)
    e_diff <- max(abs(sort(sg * cp$eigen, decreasing = TRUE) - ce)) / max(abs(ce))
    # 軌跡(千野は対象ごとに最大の長さで割って印字・平面ごとに符号の不定)
    s <- tr[tr$group == g & tr$viewpoint == vp, ]; tdiff <- 0
    for (o in unique(s$object)) { A <- matrix(cp$trajectory[o, , ], nrow = 3); A <- A / max(sqrt(colSums(Mod(A)^2)))
      for (d in 1:3) { z <- complex(real = s$real[s$object == o & s$dim == d], imaginary = s$imag[s$object == o & s$dim == d])
        tdiff <- max(tdiff, min(max(Mod(z - A[d, ])), max(Mod(z + A[d, ])))) } }
    # Lp ノルム(千野は小数2桁)
    l <- lp[lp$group == g & lp$viewpoint == vp, ]
    mine <- mapply(function(o, a, b) cp$lpnorm[a, b, o], l$object, l$time_i, l$time_j)
    ldiff <- max(abs(mine - l$lpnorm)); lrel <- max(abs(mine - l$lpnorm) / pmax(l$lpnorm, 1))
    res <- rbind(res, data.frame(group = g, view = vp, flipped = flip, eigen_reldiff = signif(e_diff, 2), traj_maxdiff = round(tdiff, 4), lp_maxdiff = round(ldiff, 3), lp_reldiff = signif(lrel, 2), n_traj = nrow(s), n_lp = nrow(l)))
  }
  # 「各時点のエルミート行列と全体のエルミート行列の距離」の正体を探す
  d <- dst$distance[dst$group == g & dst$viewpoint == 1]; Hs <- lapply(Xs, .hermitian); Xi <- f$compromise[[1]]$matrix; w <- f$compromise[[1]]$weights
  dist_mine <- f$compromise[[1]]$distance; cat(sprintf("g%d 距離の最大差 %.4f\n", g, max(abs(dist_mine - dst$distance[dst$group == g & dst$viewpoint == 1])))); next
  cand <- rbind(frob_raw = sapply(Hs, function(H) sqrt(sum(Mod(H - Xi)^2))), frob_mean = sapply(Hs, function(H) sqrt(sum(Mod(H - Xi / sum(w))^2))))
  cat(sprintf("g%d 距離 千野=%s | Frobenius(H-Xi)=%s | Frobenius(H-Xi/sum(w))=%s\n", g, paste(d, collapse = " "), paste(round(cand[1, ], 3), collapse = " "), paste(round(cand[2, ], 3), collapse = " ")))
}
print(res, row.names = FALSE)
