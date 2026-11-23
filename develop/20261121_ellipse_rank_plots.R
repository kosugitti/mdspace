# 円モデルと楕円モデルの布置を描く(職業・国の印象2ソース・車の乗り換え) --------------------------------
# MCMC のサンプルごとに，最も対数事後密度の高いサンプルへ直交プロクルステス(回転・鏡映)で合わせてから要約する。
# 楕円モデルでは座標だけを回すので，楕円の軸の向きは近似(円に近い楕円では問題にならない)

pacman::p_load(tidyverse, cmdstanr, patchwork)
source("develop/20261121_ellipse_rank_datasets.R")   # prep_data()・harshman()・mod を読み込む

datasets <- list()
datasets[["職業(表1.19)"]] <- list(s = as.matrix(read.csv("tests/testthat/fixtures/okada_imaizumi_2024_table1_19.csv",
                                                          row.names = 1, check.names = FALSE)))
x <- read.csv("tests/testthat/fixtures/okada_imaizumi_2024_table1_14.csv", check.names = FALSE)
for (b in split(x, x$source)) {
  m <- as.matrix(b[, -(1:2)])
  rownames(m) <- b$from
  datasets[[paste0("国の印象・ソース", b$source[1])]] <- list(s = m)
}
data(carswitching, package = "asymmetry")
datasets[["車の乗り換え(再尺度化)"]] <- list(s = harshman(as.matrix(carswitching))$rescaleddata)

summarise_fit <- function(pr, use_shape) {
  f <- mod$sample(data = modifyList(pr$dat, list(use_size = 1, use_shape = use_shape)),
                  init = pr$make_init, chains = 4, parallel_chains = 4, iter_warmup = 1000,
                  iter_sampling = 1000, seed = 1, refresh = 0, show_messages = FALSE,
                  output_dir = { od <- tempfile("fit_"); dir.create(od); od })
  N <- pr$dat$N
  tau <- f$draws("tau", format = "matrix")
  phi <- f$draws("phi", format = "matrix")
  lp <- as.vector(f$draws("lp__"))
  ref <- matrix(tau[which.max(lp), ], N, 2)
  # 直交プロクルステス: サンプルの布置を ref に合わせる(重心はどちらも原点)
  al <- t(apply(tau, 1, function(v) {
    Z <- matrix(v, N, 2)
    sv <- svd(t(Z) %*% ref)
    c(Z %*% sv$u %*% t(sv$v))
  }))
  tibble(obj = pr$lab, x = colMeans(al)[1:N], y = colMeans(al)[N + 1:N],
         a = apply(phi, 2, median)[1:N], b = apply(phi, 2, median)[N + 1:N])
}

ellipse_path <- function(df) {
  th <- seq(0, 2 * pi, length.out = 121)
  df %>% rowwise() %>%
    reframe(obj = obj, ex = x + a * cos(th), ey = y + b * sin(th))
}

# 表示: 円モデルでは全員の半径に同じ値を足しても mu が変わらない(水準は事前分布で決まる)ので，
# 書籍と同じく最小の半径を 0 にそろえて描く。楕円モデルは厳密な不定性ではないが，円に近いので同じ扱いにする
plot_one <- function(df, title) {
  sh <- min(pmin(df$a, df$b))
  df <- df %>% mutate(a = a - sh, b = b - sh)
  ggplot() +
    geom_path(data = ellipse_path(df), aes(ex, ey, group = obj), colour = "#3b5b8c", linewidth = 0.4) +
    geom_point(data = df, aes(x, y), size = 1.2) +
    geom_text(data = df, aes(x, y, label = obj), size = 2.4, vjust = -0.8) +
    coord_equal() +
    labs(title = title, x = "次元1", y = "次元2") +
    theme_minimal(base_family = "HiraKakuProN-W3", base_size = 9)
}

summ <- list()
for (nm in names(datasets)) {
  pr <- prep_data(datasets[[nm]]$s, "similarity")
  summ[[paste(nm, "円")]] <- summarise_fit(pr, 0)
  summ[[paste(nm, "楕円")]] <- summarise_fit(pr, 1)
}
saveRDS(summ, "develop/20261121_ellipse_rank_plots.rds")
plots <- imap(summ, function(df, key) {
  plot_one(df, sub(" 円$", ": 円(III)", sub(" 楕円$", ": 楕円(II)", key)))
})
ggsave("develop/20261121_ellipse_rank_plots.png", wrap_plots(plots, ncol = 2), width = 10, height = 19, dpi = 130)
