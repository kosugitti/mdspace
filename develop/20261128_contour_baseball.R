# 野球データ(小杉・藤原 2004)で等高線マップモデルの空間統計を試す -------------------------
# データ: ~/Dropbox/Labo/Study：バックアップ/研究：等高線/ContourMapping/邂逅/
#   野球本テスト.xls  32 人。Q1_1〜Q1_15 は類似性(対の並びは未確認)，Q2_1〜Q2_6 は好意度(5 件法)
#                    Q2 の並びは 巨人・中日・阪神・広島・ヤクルト・横浜(調査対象者 B = NO 14，C = NO 4 で確認)
#   座標など.csv      論文で使った県別の INDSCAL 共通空間(愛知 = A，兵庫 = B，広島 = C の行)
# 人の選び方は論文どおり: 大阪の 1 人と類似性に欠損のある NO 19・24 を除く(愛知 10・兵庫 9・広島 10)
# 個人のデータはリポジトリに写さず，Dropbox の元の場所から読む

pacman::p_load(tidyverse, readxl)
devtools::load_all(quiet = TRUE)
source("develop/contour_spatial.R")

base <- "~/Dropbox/Labo/Study：バックアップ/研究：等高線/ContourMapping/邂逅"
teams <- c("巨人", "中日", "阪神", "広島", "ヤクルト", "横浜")

# 評定 ------------------------------------------------------------------------
raw <- readxl::read_excel(file.path(base, "野球本テスト.xls")) %>%
  filter(live %in% c("愛知", "兵庫", "広島"), !NO %in% c(19, 24))
rat <- raw %>% select(starts_with("Q2_")) %>% as.matrix()
dimnames(rat) <- list(as.character(raw$NO), teams)

# 県別の布置(論文の図1・3・5 の布置) ----------------------------------------------------
coord <- read.csv(file.path(base, "座標など.csv"), header = FALSE, fileEncoding = "CP932")
conf_of <- function(rows) {
  m <- as.matrix(coord[rows, 2:3])
  dimnames(m) <- list(coord[rows, 1], NULL)
  return(m[teams, ])
}
confs <- list(愛知 = conf_of(1:6), 兵庫 = conf_of(8:13), 広島 = conf_of(15:20))   # 空行は読み飛ばされる

# 県ごとに contour_mds → Moran の I(W 2 種)・SAR(W 2 種) -----------------------------------
fits <- imap(confs, function(cf, pref) {
  contour_mds(cf, rat[raw$live == pref, , drop = FALSE], center = 3)
})

res <- imap_dfr(fits, function(fit, pref) {
  m_raw <- contour_moran(fit, rowstd = FALSE)
  m_std <- contour_moran(fit, rowstd = TRUE)
  s_raw <- contour_sar(fit, rowstd = FALSE)
  s_std <- contour_sar(fit, rowstd = TRUE)
  tibble(pref = pref, person = m_raw$person, n = m_raw$n,
         I_raw = m_raw$I, p_raw = m_raw$p_upper, I_std = m_std$I, p_std = m_std$p_upper,
         rho_raw = s_raw$rho, LR_raw = s_raw$LR, rho_std = s_std$rho, LR_std = s_std$LR)
})

# 県ごとのまとめ
res_pref <- res %>%
  group_by(pref) %>%
  summarise(n_person = n(), I_raw = mean(I_raw, na.rm = TRUE), I_std = mean(I_std, na.rm = TRUE),
            n_p05_raw = sum(p_raw < 0.05, na.rm = TRUE), n_p05_std = sum(p_std < 0.05, na.rm = TRUE),
            rho_std_med = median(rho_std, na.rm = TRUE), n_LR384 = sum(LR_std > 3.84, na.rm = TRUE))

# 段階1.5: 一人ずつの rho の偏り(対象 6 個・論文の 3 つの布置・真の rho を変える) -----------------
sim <- expand_grid(pref = names(confs), rho = c(-0.5, 0, 0.5, 0.8), rowstd = c(FALSE, TRUE),
                   rep = 1:500) %>%
  mutate(est = pmap_dbl(list(pref, rho, rowstd, rep), function(pref, rho, rowstd, rep) {
    set.seed(rep)
    W <- contour_W(confs[[pref]], rowstd)
    hi <- 1 / max(Re(eigen(W, only.values = TRUE)$values))
    if (rho >= hi) return(NA_real_)
    v <- sar_simulate(W, rho)
    fit <- contour_mds(confs[[pref]], v)
    contour_sar(fit, rowstd = rowstd)$rho
  }))
sim_sum <- sim %>%
  group_by(pref, rowstd, rho) %>%
  summarise(mean = mean(est, na.rm = TRUE), median = median(est, na.rm = TRUE),
            sd = sd(est, na.rm = TRUE), .groups = "drop")
