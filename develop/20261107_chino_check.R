# 千野先生の GSTATIS 出力(1999-12-30 のメール)との照合 2026-11-07 ----------------------
# 出力の原本: Labo/Study：バックアップ/修士論文/Gstatis/GSTATISのoutput.txt
# CSV にしたもの: develop/references/chino1999_gstatis/(git 管理外・ゼミ生の評定データを含む)
# 入力の 19 行列は socioKG(S/F/T/N ゼミ)と完全に一致する(グループ1=S・2=F・3=T・4=N)
# 結果: HFM・主軸回転・GSTATIS の第1段と第2段は一致。第3段の軌跡の式と欠測の埋め方は未特定

suppressMessages(devtools::load_all(quiet = TRUE))

# 1. 欠測のない行列(F ゼミ May・N ゼミ April)で HFM を照合する ---------------------------
inp <- read.csv("develop/references/chino1999_gstatis/chino1999_input.csv")
evc <- read.csv("develop/references/chino1999_gstatis/chino1999_eigenvalues.csv")
for (cs in list(c(2, "May"), c(4, "April"))) {
  g <- as.integer(cs[1]); tm <- cs[2]
  d <- inp[inp$group == g & inp$time == tm, ]; n <- max(d$row)
  M <- matrix(NA, n, n); M[cbind(d$row, d$col)] <- d$value
  f <- hfm_mds(M, ndim = 3, rotate = "none")
  ce <- evc[evc$group == g & evc$stage == tm, ]
  max(abs(sort(f$eigen, decreasing = TRUE) - ce$eigenvalue))   # 1e-5 程度になる
}

dir <- "develop/references/chino1999_gstatis/"
# 2. 千野の印字した全固有分解から各時点の H を組み直し，GSTATIS の第1段・第2段を照合する -------
ev <- read.csv(paste0(dir, "chino1999_eigenvalues.csv"))
vec <- read.csv(paste0(dir, "chino1999_eigenvectors.csv"))
rv <- read.csv(paste0(dir, "chino1999_rv.csv"))
tr <- read.csv(paste0(dir, "chino1999_trajectory.csv"))
# 千野の印字した固有分解から H を組み直す
rebuildH <- function(g, stage, vp = NA) {
  v <- vec[vec$group == g & vec$stage == stage & vec$type == "raw" & (is.na(vp) | vec$viewpoint %in% vp), ]
  ranks <- sort(unique(v$rank)); n <- max(v$object)
  H <- matrix(0i, n, n)
  for (r in ranks) { u <- v[v$rank == r, ]; z <- complex(real = u$real, imaginary = u$imag); H <- H + u$eigenvalue[1] * z %*% Conj(t(z)) }
  return(H)
}
for (g in 1:4) {
  times <- unique(ev$stage[ev$group == g & !(ev$stage %in% c("RV", "compromise"))])
  Xs <- lapply(times, function(tm) { H <- rebuildH(g, tm); X <- Re(H) + Im(H); dimnames(X) <- list(1:nrow(X), 1:nrow(X)); X })
  names(Xs) <- times
  f <- gstatis_mds(Xs, ndim = 3, ncomp = 3, rotate = "std")
  # 第1段
  crv <- rv[rv$group == g, ]; Rc <- matrix(crv$rv, length(times), byrow = TRUE)
  cev <- ev[ev$group == g & ev$stage == "RV", ]
  cat(sprintf("\n== group %d (%d 人 x %d 時点)\n", g, nrow(Xs[[1]]), length(times)))
  cat(sprintf("  RV 行列(千野は小数2桁): 最大差 %.4f\n", max(abs(round(f$rv, 2) - Rc))))
  cat(sprintf("  RV の固有値: 最大差 %.5f  (千野 %s)\n", max(abs(f$rv_eigen$values - cev$eigenvalue)), paste(cev$eigenvalue, collapse = ", ")))
  cvec <- vec[vec$group == g & vec$stage == "RV" & vec$rank == 1, ]
  cat(sprintf("  RV 第1固有ベクトル: 最大差 %.5f\n", max(abs(abs(f$rv_eigen$vectors[, 1]) - abs(cvec$real)))))
  # 第2段: 千野は重みに生の固有ベクトルを使う → こちらの固有値 x sum(v) と比べる
  for (k in 1:3) {
    ce <- sort(ev$eigenvalue[ev$group == g & ev$stage == "compromise" & ev$viewpoint == k], decreasing = TRUE)
    mine <- sort(f$compromise[[k]]$eigen, decreasing = TRUE)
    s <- if (k == 1) sum(f$rv_eigen$vectors[, 1]) else 1
    # 第2・第3視点は固有ベクトルの符号が不定 → 符号を反転したものとも比べる
    d1 <- max(abs(mine * s - ce)); d2 <- max(abs(sort(-mine * s, decreasing = TRUE) - ce))
    cat(sprintf("  視点%d 妥協の固有値: 最大差 %.5f (千野の第1固有値 %.5f)\n", k, min(d1, d2), ce[1]))
  }
}

