# パラメタリカバリ(条件を厳しくしていき，次元を選べなくなるところを探す) --------------------------
# 20261122_ellipse_rank_recovery_oi1987.R(原典の条件)では，周辺尤度が 12 データすべてで真の次元 2 を選んだ。
# 11/21 の厳しい条件(SD の比で E = 0.47〜0.77・N = 10)では選べなかったので，その間を埋める。
# 誤差 E(SD の比)= 0.5, 0.7, 1.0 と対象の数 N = 12, 10 を振り，R = 0.3 に固定。各 2 回・1〜3 次元
# 周辺尤度・順位の一致率・sigma と，circle_mds() の stress の肘のどちらが先に崩れるかを比べる

pacman::p_load(tidyverse, cmdstanr, posterior, bridgesampling)
source("develop/20261121_ellipse_rank_datasets.R")   # prep_data()・mod
devtools::load_all(quiet = TRUE)                     # circle_mds()

# 原典の条件のスクリプトから関数(simulate_oi・fit_bayes・fit_stress)だけを読み込む ------------------
local({
  ex <- parse("develop/20261122_ellipse_rank_recovery_oi1987.R")
  is_fun <- vapply(ex, function(e) {
    is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) &&
      identical(e[[3]][[1]], as.name("function"))
  }, logical(1))
  for (e in ex[is_fun]) eval(e, envir = globalenv())
})

# 条件を振って当てはめる ---------------------------------------------------------------------------
grid <- expand_grid(N = c(12, 10), E = c(0.5, 0.7, 1.0), rep = 1:2) %>%
  mutate(R = 0.3, seed = 2000 + row_number())
res <- grid %>%
  mutate(out = pmap(list(N, R, E, seed), function(N, R, E, seed) {
    sim <- simulate_oi(N = N, R = R, E = E, seed = seed)
    map_dfr(1:3, function(P) bind_cols(tibble(P = P), fit_bayes(sim, P), fit_stress(sim, P)))
  })) %>%
  unnest(out)
saveRDS(res, "develop/20261122_ellipse_rank_recovery_harsh.rds")
