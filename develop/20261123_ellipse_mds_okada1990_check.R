# ellipse_mds() を Okada (1990) の応用例(12 車種セグメントの乗り換え)と照合する -----------------------------
# 原典(develop/references/Okada_1990_generalization_scan.pdf・p.132):
#   Harshman et al. (1982) の 16 車種からスペシャルティの 4 つ(SMALL SPECIALTY DOMESTIC・CAPTIVE IMPORTS・IMPORTS，
#   MIDSIZE SPECIALTY)を除いた 12 車種。行和 + 列和(対角要素を含む)が全車種で等しくなるよう再尺度化し，
#   楕円モデルを 5〜1 次元で当てた S(stress 式 2)は 0.131, 0.155, 0.197, 0.240, 0.475
#   式 (12) の対称部と歪対称部の重みは 2 次元の解で 0.805 と 0.195
# 原典は半軸をそろえた円から始める最急降下法で，形の縛りは無い。ellipse_mds() は BFGS で複数の初期値から解くので，
# 同じ極小に止まるとは限らない(S が原典以下になるかを見る)

pacman::p_load(tidyverse)
devtools::load_all(quiet = TRUE)
source_harshman <- function() {
  for (e in parse("develop/20261121_ellipse_rank_datasets.R")) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("harshman"))) {
      eval(e, envir = globalenv())
    }
  }
}
source_harshman()

# データ: 16 車種から 4 つを除き，12 車種で再尺度化 ---------------------------------------------------
data(carswitching, package = "asymmetry")
car12 <- as.matrix(carswitching)
car12 <- car12[!(rownames(car12) %in% c("SMAD", "SMAC", "SMAI", "MIDS")),
               !(colnames(car12) %in% c("SMAD", "SMAC", "SMAI", "MIDS"))]
resc <- list("対角を含めて再尺度化(原典)" = harshman(car12, usediag = 1)$rescaleddata,
             "対角を除いて再尺度化" = harshman(car12, usediag = 0)$rescaleddata)

# 楕円モデル(II)と円モデル(III)を 5〜1 次元で当てる --------------------------------------------------
okada <- c(`5` = 0.131, `4` = 0.155, `3` = 0.197, `2` = 0.240, `1` = 0.475)
res <- expand_grid(rescale = names(resc), model = c("II", "III"), ndim = 5:1) %>%
  mutate(fit = pmap(list(rescale, model, ndim), function(rescale, model, ndim) {
    suppressWarnings(ellipse_mds(resc[[rescale]], ndim = ndim, model = model, type = "similarity",
                                 ndim_max = 5, check_orientation = FALSE))
  })) %>%
  mutate(stress = map_dbl(fit, "stress"),
         okada_1990 = okada[as.character(ndim)],
         # 半軸の範囲(円モデルは半軸を持たないので NA)
         semiaxis_min = map_dbl(fit, function(f) if (is.null(f$semiaxis)) NA_real_ else min(f$semiaxis)),
         semiaxis_max = map_dbl(fit, function(f) if (is.null(f$semiaxis)) NA_real_ else max(f$semiaxis)))
saveRDS(res, "develop/20261123_ellipse_mds_okada1990_check.rds")
