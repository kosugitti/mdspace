# 二相三元の円モデル(circle_mds にリストを渡す)の動作確認 -----------------------------------
# 表1.14(アジア10か国・2ソース)を解き，nASYMS と並べ，ソースごとの楕円を描く
# mdspace のプロジェクトのフォルダを作業ディレクトリにして，上から1節ずつ実行する
pacman::p_load(tidyverse, pkgload)

pkgload::load_all(".")

# データの読み込み(ソースごとの 10 x 10 の類似度行列のリスト) ------------------------------------
raw <- utils::read.csv("tests/testthat/fixtures/okada_imaizumi_2024_table1_14.csv",
                       check.names = FALSE)
nations <- purrr::map(1:2, function(i) {
  m <- raw %>% dplyr::filter(source == i) %>% dplyr::select(-source, -from) %>% as.matrix()
  rownames(m) <- colnames(m)
  m
})
names(nations) <- c("Source1", "Source2")

# 2次元の解(5次元から下ろす・合理的初期値＋乱数10本) -----------------------------------------
fit <- circle_mds(nations, ndim = 2, ndim_max = 5)
fit
summary(fit)

# S を nASYMS と並べる(5〜2次元は合理的初期値の経路) ---------------------------------------------
s_table <- tibble::tibble(
  ndim = 5:2,
  circle_mds = round(fit$path, 4),
  nASYMS = c(0.400, 0.431, 0.473, 0.540)
)
s_table

# 重みを nASYMS と並べる(u は全体の大きさの決め方が違うので縦横比で比べる) -----------------------------
ns_wt <- utils::read.csv("tests/testthat/fixtures/nasyms_twomode_nations_weights.csv") %>%
  dplyr::filter(ndim == 2)
weights_table <- tibble::tibble(
  source = fit$weights$source,
  w = round(fit$weights$w, 3), w_nASYMS = ns_wt$w,
  aspect = round(fit$aspect, 3),
  aspect_nASYMS = round(pmax(ns_wt$u1, ns_wt$u2) / pmin(ns_wt$u1, ns_wt$u2), 3)
)
weights_table

# 共通の布置(点と円) ---------------------------------------------------------------------
plot(fit)

# ソースごとの布置: 点を w_i 倍し，対象 j を半軸 u_it r_j の楕円で描く ------------------------------
draw_source <- function(fit, i, lim) {
  w <- fit$weights$w[i]
  u <- as.numeric(fit$weights[i, c("u1", "u2")])
  xy <- fit$conf * w
  # 楕円の半軸: v の定義(v_jki = w_i d_jk / sqrt(sum ((x_jt - x_kt) / u_it)^2))に合わせ，
  # ソース i の空間で対象 j の楕円は半軸 w_i u_it r_j
  graphics::plot(xy, asp = 1, xlim = lim, ylim = lim, pch = 16,
                 xlab = "Dimension 1", ylab = "Dimension 2",
                 main = sprintf("%s (w = %.3f, u = %.3f, %.3f)",
                                fit$weights$source[i], w, u[1], u[2]))
  graphics::abline(h = 0, v = 0, col = "grey80", lty = 3)
  th <- seq(0, 2 * pi, length.out = 181)
  for (j in seq_len(nrow(xy))) {
    a <- w * u * fit$radius[j]
    graphics::lines(xy[j, 1] + a[1] * cos(th), xy[j, 2] + a[2] * sin(th), col = "grey40")
  }
  graphics::text(xy, labels = rownames(xy), pos = 3, cex = 0.7)
}
# 枠は両ソースの楕円がすべて収まる大きさにそろえる
ext <- purrr::map_dbl(seq_len(nrow(fit$weights)), function(i) {
  w <- fit$weights$w[i]
  u <- as.numeric(fit$weights[i, c("u1", "u2")])
  max(abs(fit$conf * w) + outer(fit$radius, w * u))
})
lim <- c(-1, 1) * max(ext) * 1.05
graphics::par(mfrow = c(1, 2))
draw_source(fit, 1, lim)
draw_source(fit, 2, lim)
graphics::par(mfrow = c(1, 1))

# 乱数の初期値それぞれの S(局所解のようす) ------------------------------------------------------
sort(round(fit$starts, 4))

