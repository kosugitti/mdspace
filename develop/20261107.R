# HFM(hfm_mds)の動作確認 2026-11-07 ----------------------------------------------
# 先生が手で動かして確かめるためのスクリプト。上から順に実行する。

rm(list = ls())
pacman::p_load(tidyverse, devtools, asymmetry)
devtools::load_all()


# 1. 2019 年の教材の 4x4 の例 ---------------------------------------------------------
# jscp_lesson2.Rmd と同じ行列
x4 <- matrix(c(0, 1, 1, 7,
               1, 0, 1, 7,
               7, 7, 0, 1,
               1, 1, 7, 0), 4, byrow = TRUE,
             dimnames = list(LETTERS[1:4], LETTERS[1:4]))

# 対称部と歪対称部
dec4 <- asym_decomp(x4)
dec4
dec4$S
dec4$K
plot(dec4)

# HFM(既定は標準型の回転)
fit4 <- hfm_mds(x4, ndim = 2)
fit4
summary(fit4)
fit4$coord

# 複素平面を1次元ずつ描く
par(mfrow = c(1, 2))
plot(fit4, dim = 1)
plot(fit4, dim = 2)
par(mfrow = c(1, 1))


# 2. 2019 年の自作関数と固有値が一致するか ---------------------------------------------
hfm_2019 <- function(data) {
  s <- (data + t(data)) / 2
  sk <- (data - t(data)) / 2
  Her <- matrix(complex(re = s, im = sk), ncol = ncol(s))
  eval <- eigen(Her)$values
  eval <- eval[order(abs(eval), decreasing = TRUE)]
  return(eval)
}
tibble(old = Re(hfm_2019(x4)), new = hfm_mds(x4, ndim = 4)$eigen)


# 3. 回転の3通りを見比べる -------------------------------------------------------------
# std     : 主軸を実軸へ(実部の二乗和が最大)
# varimax : Re(z^2) の分散が最大(点が実軸か虚軸に寄る)
# none    : eigen() が返したまま
fits <- list(std = hfm_mds(x4, ndim = 2, rotate = "std"),
             varimax = hfm_mds(x4, ndim = 2, rotate = "varimax"),
             none = hfm_mds(x4, ndim = 2, rotate = "none"))
par(mfrow = c(1, 3))
plot(fits$std, dim = 1)
plot(fits$varimax, dim = 1)
plot(fits$none, dim = 1)
par(mfrow = c(1, 1))

# 回転しても H の再現は変わらない(差の最大値がほぼ 0 になるはず)
recon <- function(fit) {
  Z <- fit$coord
  return(Z %*% diag(sign(fit$eigen[seq_len(ncol(Z))]), ncol(Z)) %*% Conj(t(Z)))
}
max(Mod(recon(fits$std) - recon(fits$none)))
max(Mod(recon(fits$varimax) - recon(fits$none)))


# 4. 2001 年の VB 版(rotateSTD)の向きの確認 -----------------------------------------
# VB は z * exp(+i th) で回していた。主軸を実軸へそろえるには exp(-i th) が要るはず
z <- fits$none$coord[, 1]
A <- sum(Re(z)^2 - Im(z)^2)
B <- 2 * sum(Re(z) * Im(z))
th <- (atan(B / A) + ifelse(A < 0, pi, 0)) / 2
tibble(rotation = c("VB (+th)", "逆向き (-th)", "mdspace std"),
       sum_re2 = c(sum(Re(z * exp(1i * th))^2),
                   sum(Re(z * exp(-1i * th))^2),
                   sum(Re(fits$std$coord[, 1])^2)))


# 5. 実データ1: asymmetry パッケージの学生の移動 -----------------------------------
data(studentmigration, package = "asymmetry")
mig <- as_proximity(studentmigration, type = "similarity")
asym_decomp(mig)
plot(asym_decomp(mig))

fit_mig <- hfm_mds(mig, ndim = 3)
fit_mig
plot(fit_mig, dim = 1)
plot(fit_mig, dim = 2)


# 6. 実データ2: socioKG(KAKEN24730510 のゼミの関係評定・N ゼミ 10 人の T1) ----------
# 2014 年の gstatis の例と同じく，欠測は 4 で埋める
load("~/Library/CloudStorage/Dropbox/Labo/Archive：研究会資料/CognitiveMiser/KAKEN24730510/data/socioKG.rda")
nT1 <- socioKG$Nsemi$T1
nT1[is.na(nT1)] <- 4
dimnames(nT1) <- list(paste0("P", 1:10), paste0("P", 1:10))

fit_soc <- hfm_mds(nT1, ndim = 2)
fit_soc
par(mfrow = c(1, 2))
plot(fit_soc, dim = 1)
plot(fit_soc, dim = 2)
par(mfrow = c(1, 1))



# 7. GSTATIS: socioKG の N ゼミ 10 人 × 5 時点 -------------------------------------
# 欠測は 4 で埋める(2014 年の例と同じ)
socio_list <- socioKG$Nsemi[paste0("T", 1:5)] %>%
  map(function(m) {
    m[is.na(m)] <- 4
    dimnames(m) <- list(paste0("P", 1:10), paste0("P", 1:10))
    return(m)
  })

fit_gst <- gstatis_mds(socio_list, ndim = 2)
fit_gst

# 第1段: 時点どうしの RV 係数と，時点の配置
round(fit_gst$rv, 3)
plot(fit_gst, type = "interstructure")

# 第2段: 妥協(5 時点の加重平均)の HFM
plot(fit_gst, type = "compromise", dim = 1)
plot(fit_gst, type = "compromise", dim = 2)

# 第3段: 各時点を妥協の平面へ射影した軌跡(黒丸が妥協の位置，T5 が最後の時点)
par(mfrow = c(1, 2))
plot(fit_gst, type = "trajectory", dim = 1)
plot(fit_gst, type = "trajectory", dim = 2)
par(mfrow = c(1, 1))

# 視点を2つ作る(RV の第2成分を重みにした妥協)
fit_gst2 <- gstatis_mds(socio_list, ndim = 2, ncomp = 2)
fit_gst2$weights
plot(fit_gst2, type = "compromise", dim = 1, view = 2)

# 時点ごとに行列の大きさをそろえる場合(STATIS の通常の正規化)
fit_gst_norm <- gstatis_mds(socio_list, ndim = 2, scale = "norm")
fit_gst_norm$weights
