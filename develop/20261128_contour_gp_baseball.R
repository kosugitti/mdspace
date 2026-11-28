# ガウス過程版の等高線マップモデル(develop/contour_gp.stan)を野球データで試す(2026-11-28)
# データと人の選び方は 20261128_contour_baseball.R と同じ。4 本の連鎖を 1 コアで順に走らせる(23 秒)

suppressMessages({pacman::p_load(tidyverse, readxl, cmdstanr, posterior); devtools::load_all(quiet = TRUE)})
base <- "~/Dropbox/Labo/Study：バックアップ/研究：等高線/ContourMapping/邂逅"
teams <- c("巨人", "中日", "阪神", "広島", "ヤクルト", "横浜")
raw <- readxl::read_excel(file.path(base, "野球本テスト.xls")) %>% filter(live %in% c("愛知","兵庫","広島"), !NO %in% c(19,24))
rat <- as.matrix(raw[, paste0("Q2_",1:6)]) - 3
coord <- read.csv(file.path(base, "座標など.csv"), header = FALSE, fileEncoding = "CP932")
conf_of <- function(rows) { m <- as.matrix(coord[rows, 2:3]); rownames(m) <- coord[rows,1]; m[teams, ] }
prefs <- c("愛知","兵庫","広島"); confs <- list(conf_of(1:6), conf_of(8:13), conf_of(15:20))
long <- data.frame(i = rep(seq_len(nrow(rat)), each = 6), obj = rep(1:6, nrow(rat)), v = as.vector(t(rat))) %>% filter(!is.na(v))
st <- long %>% mutate(r = row_number()) %>% group_by(i) %>% summarise(start = min(r), len = n())
dat <- list(G = 3, J = 6, X = confs, N = nrow(rat), grp = match(raw$live, prefs), L = nrow(long),
            v = long$v, obj = long$obj, start = st$start, len = st$len)
mod <- cmdstan_model("develop/contour_gp.stan")
t0 <- Sys.time()
fit <- mod$sample(data = dat, chains = 4, parallel_chains = 1, iter_warmup = 1000, iter_sampling = 1000,
                  seed = 1, refresh = 0, adapt_delta = 0.95, show_messages = FALSE)
cat("秒:", round(as.numeric(difftime(Sys.time(), t0, units="secs")),1), "\n")
s <- fit$summary(c("ell","eta0","m_eta","s_eta","sig","alpha0","s_alpha"))
print(s %>% select(variable, mean, q5, q95, rhat, ess_bulk), digits = 3)
print(fit$diagnostic_summary())
om <- fit$summary("omega") %>% mutate(NO = raw$NO, pref = raw$live, concern = raw$concern)
print(om %>% select(NO, pref, concern, mean, q5, q95), n = 30, digits = 2)
f <- fit$summary("f") ; print(f %>% select(variable, mean, q5, q95), n = 18, digits = 2)
