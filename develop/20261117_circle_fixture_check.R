# 円モデル(circle_mds)の下準備: 表1.19と nASYMS の2次元解の整合を確かめる -----------------
pacman::p_load(tidyverse)

# データと nASYMS の解の読み込み ------------------------------------------------
X <- utils::read.csv("tests/testthat/fixtures/okada_imaizumi_2024_table1_19.csv",
                     row.names = 1, check.names = FALSE) %>% as.matrix()
conf <- utils::read.csv("tests/testthat/fixtures/nasyms_radius_occupation_2d.csv")
Z <- as.matrix(conf[, c("dim1", "dim2")])
r <- conf$radius

# m_jk = d_jk - r_j + r_k を作る ---------------------------------------------------
D <- as.matrix(stats::dist(Z))
M <- D - outer(r, r, "-")
off <- row(X) != col(X)
dat <- tibble::tibble(s = X[off], m = M[off], j = row(X)[off], k = col(X)[off])

# マニュアル Figure 10 の先頭3行(m = 0.259, 0.06, 0.044)と照合 --------------------------
check_rows <- dat %>% dplyr::arrange(dplyr::desc(s)) %>% dplyr::slice(1:3)

# 類似度なので m は s の単調減少で当てる(同順位は一次の扱い=順位の付け替え可を省略) ------------
ord <- order(-dat$s)
mhat <- numeric(nrow(dat))
mhat[ord] <- stats::isoreg(dat$m[ord])$yf
S1 <- sqrt(sum((dat$m - mhat)^2) / sum(dat$m^2))
S2 <- sqrt(sum((dat$m - mhat)^2) / sum((dat$m - mean(dat$m))^2))
res <- list(check_rows = check_rows, mhat_top3 = mhat[ord][1:3], S1 = S1, S2 = S2)
res
