# 足場(点1〜3)の選び方で事後分布が変わらないかを確かめる(Stan, MCMC) -------------------------
# 事前分布を重心からのずれに置き，回転の縛りの補正 log(a) を入れた vonmises_mds.stan で，
# 千野データを2通りの足場で回し，重心を原点にした布置の点間距離と mu の事後平均を比べる

pacman::p_load(tidyverse, cmdstanr, posterior)
source("develop/20261118_vonmises_stan.R", echo = FALSE)   # mod・chino・fit の関数(全体が走る)

fit_with_order <- function(G, ord, seed = 1) {
  Gs <- G[ord, ord]
  diag(Gs) <- NA
  sc <- sqrt(mean(Gs^2, na.rm = TRUE))
  use <- !is.na(Gs)
  dat <- list(N = nrow(G), L = sum(use), row = row(Gs)[use], col = col(Gs)[use], g = Gs[use] / sc,
              omega_fixed = 1, kappa_max = 6.5, prior_x_sd = 1, prior_c_a = 2, prior_c_b = 4)
  opts <- map(1:10, ~ tryCatch(mod$optimize(data = dat, seed = .x, refresh = 0, iter = 5000,
                                            show_messages = FALSE, show_exceptions = FALSE),
                               error = function(e) NULL)) %>% compact()
  best <- opts[[which.max(map_dbl(opts, ~ .x$lp()))]]
  init <- map(1:4, ~ list(Lx = matrix(unname(best$mle("Lx")), nrow(G) - 1, 2),
                          m = map(unname(best$mle("mu")), function(a) c(cos(a), sin(a))),
                          kappa = pmin(unname(best$mle("kappa")), 6.4), sigma = unname(best$mle("sigma"))))
  fit <- mod$sample(data = dat, seed = seed, chains = 4, parallel_chains = 4, iter_warmup = 1000,
                    iter_sampling = 3000, refresh = 0, show_messages = FALSE, init = init)
  # 重心を原点にした布置から点間距離(回転によらない)を事後平均で，元の並びに戻す
  Xc <- fit$draws("Xc", format = "draws_matrix")
  n <- nrow(G)
  back <- order(ord)
  D <- apply(Xc, 1, function(v) as.vector(as.matrix(dist(matrix(v, n, 2)[back, ]))))
  list(D = matrix(rowMeans(D), n, n, dimnames = dimnames(G)),
       kappa = colMeans(fit$draws("kappa", format = "draws_matrix"))[back],
       rhat = max(fit$summary(c("Xc", "kappa"))$rhat, na.rm = TRUE))
}

ordA <- choose_anchor(chino)
# 足場 B: 足場 A の3点を除いた7人から，同じ基準(離れていて一直線に並ばない)で選び直す
rest <- setdiff(seq_len(10), ordA[1:3])
subB <- choose_anchor(chino[rest, rest])
ordB <- c(rest[subB[1:3]], setdiff(seq_len(10), rest[subB[1:3]]))
fA <- fit_with_order(chino, ordA)
fB <- fit_with_order(chino, ordB)
tibble(anchor_A = paste(rownames(chino)[ordA[1:3]], collapse = ""),
       anchor_B = paste(rownames(chino)[ordB[1:3]], collapse = ""),
       maxdiff_dist = max(abs(fA$D - fB$D)), mean_dist = mean(fA$D[upper.tri(fA$D)]),
       maxdiff_kappa = max(abs(fA$kappa - fB$kappa)), rhat_A = fA$rhat, rhat_B = fB$rhat)
