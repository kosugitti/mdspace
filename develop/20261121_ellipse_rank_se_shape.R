# 形の散らばり s_e の事後分布の形を見る(MCMC と最適化・Pathfinder の食い違いの確かめ) ---------------
pacman::p_load(tidyverse, cmdstanr, patchwork)

fit <- readRDS("develop/20261121_ellipse_rank_fit.rds")
vbr <- readRDS("develop/20261121_ellipse_rank_vb.rds")
dr <- fit$draws(c("s_e", "s_c", "sig", "lp__"), format = "df") %>% as_tibble()
map_se <- vbr$tab %>% filter(par == "s_e") %>% pull(Optimize)
pf_se <- vbr$tab %>% filter(par == "s_e") %>% pull(Pathfinder)

ink <- "#3b5b8c"
p1 <- ggplot(dr, aes(s_e)) +
  geom_histogram(bins = 60, fill = ink, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = map_se, linetype = 2) +
  geom_vline(xintercept = pf_se, linetype = 3) +
  annotate("text", x = map_se, y = Inf, label = "最適化", vjust = 1.5, hjust = -0.1, size = 3) +
  annotate("text", x = pf_se, y = Inf, label = "Pathfinder", vjust = 3, hjust = -0.1, size = 3) +
  labs(x = "s_e(形の散らばり)", y = "度数", title = "(a) s_e の事後分布(MCMC 4,000)") +
  theme_minimal(base_family = "HiraKakuProN-W3")
p2 <- ggplot(dr, aes(s_e, s_c)) +
  geom_point(alpha = 0.25, size = 0.8, colour = ink) +
  labs(x = "s_e(形の散らばり)", y = "s_c(大きさの散らばり)", title = "(b) s_e と s_c") +
  theme_minimal(base_family = "HiraKakuProN-W3")
p3 <- ggplot(dr, aes(s_e, lp__)) +
  geom_point(alpha = 0.25, size = 0.8, colour = ink) +
  labs(x = "s_e(形の散らばり)", y = "lp__(対数事後密度・ヤコビアンあり)", title = "(c) s_e と対数事後密度") +
  theme_minimal(base_family = "HiraKakuProN-W3")
ggsave("develop/20261121_ellipse_rank_se_shape.png", p1 / (p2 | p3), width = 9, height = 7, dpi = 150)

qs <- quantile(dr$s_e, c(0.05, 0.25, 0.5, 0.75, 0.95))
saveRDS(list(q = qs, p_gt_08 = mean(dr$s_e > 0.8), p_lt_01 = mean(dr$s_e < 0.1),
             by_chain = fit$draws("s_e", format = "df") %>% as_tibble() %>%
               group_by(.chain) %>% summarise(mean = mean(s_e), q95 = quantile(s_e, 0.95))),
        "develop/20261121_ellipse_rank_se_shape.rds")
