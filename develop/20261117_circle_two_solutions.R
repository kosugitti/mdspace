# 職業データ2次元: S = 0.3647 の解と 0.3716 の解(nASYMS と同じ谷)を比べる -----------------
pacman::p_load(tidyverse, pkgload)

pkgload::load_all(".")

# データと2つの解 -------------------------------------------------------------------
occ <- utils::read.csv("tests/testthat/fixtures/okada_imaizumi_2024_table1_19.csv",
                       row.names = 1, check.names = FALSE) %>%
  as.matrix()
fit_best <- circle_mds(occ, ndim = 2, ndim_max = 5)               # 乱数込みの最良
fit_rat <- circle_mds(occ, ndim = 2, ndim_max = 5, nrandom = 0)   # 合理的初期値だけ

# 最良の解を合理的初期値の解に直交回転(鏡映込み)でそろえる ------------------------------------
sv <- svd(t(fit_best$conf) %*% fit_rat$conf)
conf_best <- fit_best$conf %*% sv$u %*% t(sv$v)
dimnames(conf_best) <- dimnames(fit_rat$conf)

# 座標と半径の表 ----------------------------------------------------------------------
compare <- tibble::tibble(
  object = rownames(occ),
  D1_rational = round(fit_rat$conf[, 1], 3), D1_best = round(conf_best[, 1], 3),
  D2_rational = round(fit_rat$conf[, 2], 3), D2_best = round(conf_best[, 2], 3),
  r_rational = round(fit_rat$radius, 3), r_best = round(fit_best$radius, 3)
)
compare

# 2つの解の図を並べる -------------------------------------------------------------------
draw <- function(conf, radius, main, lim) {
  graphics::plot(conf, asp = 1, xlim = lim, ylim = lim, pch = 16,
                 xlab = "Dimension 1", ylab = "Dimension 2", main = main)
  graphics::abline(h = 0, v = 0, col = "grey80", lty = 3)
  th <- seq(0, 2 * pi, length.out = 181)
  for (i in seq_len(nrow(conf))) {
    graphics::lines(conf[i, 1] + radius[i] * cos(th),
                    conf[i, 2] + radius[i] * sin(th), col = "grey40")
  }
  graphics::text(conf, labels = rownames(conf), pos = 3, cex = 0.7)
}
lim <- range(c(fit_rat$conf - fit_rat$radius, fit_rat$conf + fit_rat$radius,
               conf_best - fit_best$radius, conf_best + fit_best$radius))
grDevices::png("develop/sim_results/20261117_circle_two_solutions.png",
               width = 1600, height = 800, res = 130)
graphics::par(mfrow = c(1, 2))
draw(fit_rat$conf, fit_rat$radius,
     sprintf("Rational start (S = %.4f)", fit_rat$stress), lim)
draw(conf_best, fit_best$radius,
     sprintf("Best of 11 starts (S = %.4f)", fit_best$stress), lim)
grDevices::dev.off()
