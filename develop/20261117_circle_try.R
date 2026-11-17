# circle_mds() の動作確認(書籍表1.19・nASYMS との比較) -------------------------------
# mdspace のプロジェクトのフォルダを作業ディレクトリにして実行する
pacman::p_load(tidyverse, pkgload)

pkgload::load_all(".")

# データの読み込み ----------------------------------------------------------------
occ <- utils::read.csv("tests/testthat/fixtures/okada_imaizumi_2024_table1_19.csv",
                       row.names = 1, check.names = FALSE) %>%
  as.matrix()
nasyms <- utils::read.csv("tests/testthat/fixtures/nasyms_radius_occupation_all.csv")

# 2次元の解(5次元から1次元ずつ下ろす・nASYMS と同じ合理的初期値だけ) ------------------------
fit2 <- circle_mds(occ, ndim = 2, ndim_max = 5, nrandom = 0)
fit2
summary(fit2)
plot(fit2)

# 5次元から1次元までの S を nASYMS と並べる --------------------------------------------
fit1 <- circle_mds(occ, ndim = 1, ndim_max = 5, nrandom = 0)
s_table <- tibble::tibble(
  ndim = 5:1,
  circle_mds = round(fit1$path, 4),
  nASYMS = c(0.183, 0.202, 0.276, 0.372, 0.499)
)
s_table

# 2次元の布置と半径を nASYMS と比べる(直交回転をそろえる) -------------------------------
ns2 <- nasyms %>% dplyr::filter(ndim == 2)
Z <- as.matrix(ns2[, c("D1", "D2")])
sv <- svd(t(fit2$conf) %*% Z)
conf_rot <- fit2$conf %*% sv$u %*% t(sv$v)
compare2 <- tibble::tibble(
  object = ns2$object,
  D1 = round(conf_rot[, 1], 3), D1_nASYMS = ns2$D1,
  D2 = round(conf_rot[, 2], 3), D2_nASYMS = ns2$D2,
  radius = round(fit2$radius, 3), radius_nASYMS = ns2$radius
)
compare2

# 既定(合理的初期値＋乱数の初期値10本)で各次元を解く ---------------------------------------
fits <- purrr::map(1:5, ~ circle_mds(occ, ndim = .x, ndim_max = 5))
s_multi <- tibble::tibble(
  ndim = 1:5,
  rational = purrr::map_dbl(fits, ~ round(.x$starts[1], 4)),
  best = purrr::map_dbl(fits, ~ round(.x$stress, 4)),
  best_start = purrr::map_int(fits, "best_start"),
  nASYMS = c(0.499, 0.372, 0.276, 0.202, 0.183)
)
s_multi
fits[[2]]
plot(fits[[2]])
# 2次元: 11本の初期値それぞれの S
sort(round(fits[[2]]$starts, 4))
