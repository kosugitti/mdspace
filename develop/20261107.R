# HFM(hfm_mds)の動作確認 2026-11-07 ----------------------------------------------
# 先生が手で動かして確かめるためのスクリプト。上から順に実行する。

rm(list = ls())
pacman::p_load(tidyverse, devtools, asymmetry)
devtools::load_all()


# 1. 2019 年の教材の 4x4 の例 ---------------------------------------------------------
# jscp_lesson2.Rmd と同じ行列
x4 <- matrix(
  c(
    0, 1, 1, 7,
    1, 0, 1, 7,
    7, 7, 0, 1,
    1, 1, 7, 0
  ), 4,
  byrow = TRUE,
  dimnames = list(LETTERS[1:4], LETTERS[1:4])
)

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
fits <- list(
  std = hfm_mds(x4, ndim = 2, rotate = "std"),
  varimax = hfm_mds(x4, ndim = 2, rotate = "varimax"),
  none = hfm_mds(x4, ndim = 2, rotate = "none")
)
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
tibble(
  rotation = c("VB (+th)", "逆向き (-th)", "mdspace std"),
  sum_re2 = c(
    sum(Re(z * exp(1i * th))^2),
    sum(Re(z * exp(-1i * th))^2),
    sum(Re(fits$std$coord[, 1])^2)
  )
)


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

# 第2段: 妥協(GRV の第1固有ベクトルを重みにした5時点の和)の HFM
plot(fit_gst, type = "compromise", dim = 1)
plot(fit_gst, type = "compromise", dim = 2)

# 第3段: 各時点の行ベクトルを妥協の平面へ射影した軌跡
# (白丸が各時点，黒丸が重みつき平均の位置，T5 が最後の時点。scale = "max" で千野版と同じ縮尺)
par(mfrow = c(1, 2))
plot(fit_gst, type = "trajectory", dim = 1)
plot(fit_gst, type = "trajectory", dim = 2)
par(mfrow = c(1, 1))
plot(fit_gst, type = "trajectory", dim = 1, objects = c("P2", "P4", "P8"), scale = "max")

# 時点間の Lp ノルム(対象 P2)と，各時点の H と妥協 Xi の距離
round(fit_gst$lpnorm[, , "P2"], 2)
round(fit_gst$compromise$C1$distance, 2)

# 視点を2つ作る(RV の第2成分を重みにした妥協)
fit_gst2 <- gstatis_mds(socio_list, ndim = 2, ncomp = 2)
fit_gst2$weights
plot(fit_gst2, type = "compromise", dim = 1, view = 2)

# 時点ごとに行列の大きさをそろえる場合(STATIS の通常の正規化)
fit_gst_norm <- gstatis_mds(socio_list, ndim = 2, scale = "norm")
fit_gst_norm$weights


# 8. 千野先生の FORTRAN 版(1999)の出力と目で見比べる ------------------------------------
# 欠測の埋め方が特定できないので，千野が印字した全固有分解から補完後の行列を組み直して使う
chino_dir <- "develop/references/chino1999_gstatis/"
ch_ev <- read_csv(paste0(chino_dir, "chino1999_eigenvalues.csv"), show_col_types = FALSE)
ch_vec <- read_csv(paste0(chino_dir, "chino1999_eigenvectors.csv"), show_col_types = FALSE)
ch_rv <- read_csv(paste0(chino_dir, "chino1999_rv.csv"), show_col_types = FALSE)
ch_tr <- read_csv(paste0(chino_dir, "chino1999_trajectory.csv"), show_col_types = FALSE)
ch_lp <- read_csv(paste0(chino_dir, "chino1999_lpnorm.csv"), show_col_types = FALSE)
ch_dist <- read_csv(paste0(chino_dir, "chino1999_distance.csv"), show_col_types = FALSE)

# 千野の固有分解 H = sum(lambda u u*) から元の行列 X = Re(H) + Im(H) を組み直す
chino_matrix <- function(g, stage) {
  v <- ch_vec %>% filter(group == g, stage == !!stage, type == "raw")
  n <- max(v$object)
  H <- matrix(0i, n, n)
  for (r in unique(v$rank)) {
    u <- v %>% filter(rank == r)
    z <- complex(real = u$real, imaginary = u$imag)
    H <- H + u$eigenvalue[1] * z %*% Conj(t(z))
  }
  X <- Re(H) + Im(H)
  dimnames(X) <- list(1:n, 1:n)
  return(X)
}

# グループを選ぶ(1=S ゼミ 21人x4, 2=F ゼミ 14人x5, 3=T ゼミ 20人x5, 4=N ゼミ 10人x5)
g <- 1
times_g <- ch_ev %>%
  filter(group == g, !stage %in% c("RV", "compromise")) %>%
  pull(stage) %>%
  unique()
xs_g <- map(times_g, ~ chino_matrix(g, .x)) %>% set_names(times_g)
fit_ch <- gstatis_mds(xs_g, ndim = 3, ncomp = 3)

# 8-1. GRV 行列(千野は小数2桁)
ch_rv %>%
  filter(group == g) %>%
  mutate(mdspace = map2_dbl(time_i, time_j, ~ fit_ch$rv[.x, .y])) %>%
  mutate(mdspace_2digits = round(mdspace, 2))

# 8-2. GRV の固有値と第1固有ベクトル(=妥協の重み)
tibble(
  rank = seq_along(times_g),
  chino = ch_ev %>% filter(group == g, stage == "RV") %>% pull(eigenvalue),
  mdspace = fit_ch$rv_eigen$values
)
tibble(
  time = times_g,
  chino = ch_vec %>% filter(group == g, stage == "RV", rank == 1) %>% pull(real),
  mdspace = fit_ch$weights[, 1]
)

# 8-3. 視点1の妥協の固有値(上位6つ)
tibble(
  chino = ch_ev %>% filter(group == g, stage == "compromise", viewpoint == 1) %>%
    pull(eigenvalue) %>% sort(decreasing = TRUE) %>% head(6),
  mdspace = sort(fit_ch$compromise$C1$eigen, decreasing = TRUE) %>% head(6)
)

# 8-4. 各時点の H と妥協 Xi の距離
ch_dist %>%
  filter(group == g, viewpoint == 1) %>%
  mutate(mdspace = fit_ch$compromise$C1$distance[time])

# 8-5. Lp ノルム: 千野(横軸)と mdspace(縦軸)。対角線に乗れば一致
lp_cmp <- ch_lp %>%
  filter(group == g, viewpoint <= 3, time_i != time_j) %>%
  mutate(mdspace = pmap_dbl(
    list(viewpoint, object, time_i, time_j),
    function(v, o, a, b) fit_ch$compromise[[v]]$lpnorm[a, b, o]
  ))
ggplot(lp_cmp, aes(x = lpnorm, y = mdspace, colour = factor(viewpoint))) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60") +
  geom_point(alpha = 0.6) +
  labs(
    x = "Chino (FORTRAN)", y = "mdspace", colour = "view",
    title = paste0("Lp norm: group ", g)
  ) +
  theme_minimal()

# 8-6. 軌跡: 千野(赤の x)と mdspace(青の o)を同じ平面に重ねる
# 千野は対象ごとに最大の長さで割って印字しているので，mdspace も同じ縮尺にする
# 視点2・3は重みの符号の決め方が違うことがある(千野版は手順まかせ)。ずれたら -X で回した解と比べる
traj_cmp <- function(fit, g, view = 1) {
  cp <- fit$compromise[[view]]
  ch <- ch_tr %>% filter(group == g, viewpoint == view)
  out <- ch %>%
    distinct(object) %>%
    pull(object) %>%
    map_dfr(function(o) {
      A <- matrix(cp$trajectory[o, , ], nrow = dim(cp$trajectory)[2])
      A <- A / max(sqrt(colSums(Mod(A)^2)))
      ch %>%
        filter(object == o) %>%
        mutate(my = A[cbind(dim, match(time, dimnames(cp$trajectory)[[3]]))])
    })
  # 平面ごとに符号(180 度)の不定があるので，千野に近い向きにそろえる
  out <- out %>%
    group_by(object, dim) %>%
    mutate(
      flip = sum(Mod(complex(real = real, imaginary = imag) + my)) <
        sum(Mod(complex(real = real, imaginary = imag) - my)),
      my = if_else(flip, -my, my)
    ) %>%
    ungroup()
  return(out)
}
tc <- traj_cmp(fit_ch, g, view = 1)
# 次元ごとに1枚ずつ描く(縦横比は 1)
for (d in sort(unique(tc$dim))) {
  td <- tc %>% filter(dim == d)
  p <- ggplot(td) +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_vline(xintercept = 0, colour = "grey85") +
    geom_path(aes(x = Re(my), y = Im(my), group = object), colour = "steelblue", alpha = 0.5) +
    geom_point(aes(x = Re(my), y = Im(my)), colour = "steelblue", shape = 1, size = 4) +
    geom_point(aes(x = real, y = imag), colour = "firebrick", shape = 4, size = 2) +
    geom_text(
      data = td %>% group_by(object) %>% slice_tail(n = 1),
      aes(x = real, y = imag, label = paste0("obj", object)), size = 3, vjust = -1
    ) +
    coord_equal() +
    labs(
      x = "Real", y = "Imaginary",
      title = paste0(
        "Trajectories: group ", g, ", view 1, dim ", d,
        " (red x = Chino, blue o = mdspace)"
      )
    ) +
    theme_minimal()
  print(p)
}

# 数値でも: 千野と mdspace の差の最大
tc %>% summarise(max_diff = max(Mod(complex(real = real, imaginary = imag) - my)))


# 9. 欠測の扱い(missing = "impute"・MINERS 法) -----------------------------------------
# 片側だけの欠け → 転置の値で埋める(その組の歪対称部は 0)
# 対角と両側の欠け → MINERS 法(行と列の平均から始め，HFM の上位 p 次元の再構成で置き換えを繰り返す)
# その時点にいない人(行も列も全部欠け) → データの不備としてエラー

# 9-1. socioKG の4つのゼミの欠測の型
sockg_raw <- map(socioKG[c("Ssemi", "Fsemi", "Tsemi", "Nsemi")], function(z) {
  ts <- names(z)[str_detect(names(z), "^T[0-9]$")]
  return(z[ts])
})
map_dfr(sockg_raw, missing_pattern, .id = "zemi")

# 9-2. N ゼミ: 欠測を 4 で埋めずにそのまま渡す
n_raw <- sockg_raw$Nsemi %>%
  map(function(m) {
    dimnames(m) <- list(paste0("P", 1:10), paste0("P", 1:10))
    return(m)
  })
missing_pattern(n_raw)

# 既定(missing = "none")ではエラーになる
try(hfm_mds(n_raw$T3, ndim = 2))

# T3 を MINERS 法で埋めて HFM
fit_n3 <- hfm_mds(n_raw$T3, ndim = 2, missing = "impute")
fit_n3
fit_n3$imputed
fit_n3$impute_niter

# 4 で埋めた場合と比べる
n3_four <- n_raw$T3
n3_four[is.na(n3_four)] <- 4
tibble(impute = fit_n3$eigen[1:5], fill4 = hfm_mds(n3_four, ndim = 2)$eigen[1:5])

# T2 は P2 がいない(行も列も全部欠け)ので，GSTATIS はエラーになる
try(gstatis_mds(n_raw, ndim = 2, missing = "impute"))

# P2 がいない T2 を外せば通る
fit_n_imp <- gstatis_mds(n_raw[c("T1", "T3", "T4", "T5")], ndim = 2, missing = "impute")
fit_n_imp
fit_n_imp$imputed
plot(fit_n_imp, type = "trajectory", dim = 1)

# 9-3. 千野版と比べる: 片側だけの欠け(F ゼミ July)は千野と一致するはず
g2_july <- ch_vec %>% filter(group == 2, stage == "July")   # 8 節の ch_vec を使う
f_july <- sockg_raw$Fsemi$T3
dimnames(f_july) <- list(1:14, 1:14)
tibble(chino = ch_ev %>% filter(group == 2, stage == "July") %>% pull(eigenvalue),
       mdspace = sort(hfm_mds(f_july, ndim = 3, missing = "impute")$eigen, decreasing = TRUE))
