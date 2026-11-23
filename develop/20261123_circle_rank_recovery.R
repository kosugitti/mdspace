# パラメタリカバリのやり直し(circle_rank.stan: 閾値を潜在変数の広がりに合わせ，回転を足場で止めた版) -----
# 11/22 のリカバリ(20261122_ellipse_rank_recovery_oi1987.R・_harsh.R)は，閾値を sig の単位に置いた版で，
# 閾値の事前分布が sig の強い事前分布として働いていた(E = 0.5 以上で 1 次元ばかり選んだ)。
# 同じ 24 データ(同じ乱数の種)に当て直し，周辺尤度がどの条件まで真の 2 次元を選べるかを見る
#   原典の条件: N = 12・R = 0.1, 0.3, 0.5・E = 0.1, 0.3・各 2 回(種 1001〜1012)
#   厳しい条件: N = 12, 10・R = 0.3・E = 0.5, 0.7, 1.0・各 2 回(種 2001〜2012)
# circle_mds() の stress は 11/22 の rds にあるので，ここでは計算しない

pacman::p_load(tidyverse, cmdstanr, posterior, bridgesampling)
mod_circle <- cmdstanr::cmdstan_model("develop/circle_rank.stan", compile_model_methods = TRUE,
                                      force_recompile = TRUE)

# 関数を読み込む(simulate_oi と，足場の版の pick_scaffold・to_scaffold・prep_circle・fit_circle) ----------
local({
  load_funs <- function(file, names) {
    for (e in parse(file)) {
      if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]])[1] %in% names) {
        eval(e, envir = globalenv())
      }
    }
  }
  load_funs("develop/20261122_ellipse_rank_recovery_oi1987.R", "simulate_oi")
  load_funs("develop/20261123_circle_rank_scaffold.R",
            c("pick_scaffold", "to_scaffold", "prep_circle", "fit_circle"))
})

# 条件 ----------------------------------------------------------------------------------------------
grid <- bind_rows(
  expand_grid(R = c(0.1, 0.3, 0.5), E = c(0.1, 0.3), rep = 1:2) %>%
    mutate(N = 12, seed = 1000 + row_number(), cond = "原典"),
  expand_grid(N = c(12, 10), E = c(0.5, 0.7, 1.0), rep = 1:2) %>%
    mutate(R = 0.3, seed = 2000 + row_number(), cond = "厳しい")
)
res <- grid %>%
  mutate(out = pmap(list(N, R, E, seed), function(N, R, E, seed) {
    sim <- simulate_oi(N = N, R = R, E = E, seed = seed)
    off <- row(sim$D) != col(sim$D)
    M <- sim$D - outer(sim$r, rep(1, N)) + outer(rep(1, N), sim$r)
    map_dfr(1:3, function(P) {
      bind_cols(tibble(P = P, sigma_true = E * sd(M[off])), fit_circle(sim, P))
    })
  })) %>%
  unnest(out)
saveRDS(res, "develop/20261123_circle_rank_recovery.rds")
