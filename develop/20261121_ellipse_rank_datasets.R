# 楕円モデルの順位尤度版を他のデータに当てる: 4つの下位モデル(楕円・円・等面積の楕円・対称)を ----------
# MCMC で推定し，周辺尤度(ブリッジサンプリング)で比べる。職業データと同じ手順を関数にした
# 使い方: Rscript develop/20261121_ellipse_rank_datasets.R nations

pacman::p_load(tidyverse, cmdstanr, posterior, bridgesampling)

# データを Stan に渡す形にする(類似度は符号を反転して非類似度の順位に) -------------------------------
prep_data <- function(s, type = c("similarity", "dissimilarity"), P = 2) {
  type <- match.arg(type)
  s <- as.matrix(s)
  diag(s) <- NA
  N <- nrow(s)
  use <- row(s) != col(s) & !is.na(s)
  d <- if (type == "similarity") -s[use] else s[use]
  vals <- sort(unique(d))
  # 基準の対象: 対称部の Torgerson 解で，軸ごとに重心から最も遠い対象
  D0 <- if (type == "similarity") max(s, na.rm = TRUE) - s else s
  D0[is.na(D0)] <- t(D0)[is.na(D0)]
  D0[is.na(D0)] <- mean(D0, na.rm = TRUE)
  diag(D0) <- 0
  Z0 <- stats::cmdscale((D0 + t(D0)) / 2, k = P)
  anchor <- apply(abs(sweep(Z0, 2, colMeans(Z0))), 2, which.max)
  Z0 <- sweep(Z0, 2, sign(Z0[cbind(anchor, seq_len(P))]), `*`)
  dat <- list(N = N, P = P, L = sum(use), row = row(s)[use], col = col(s)[use],
              K = length(vals), y = match(d, vals), anchor = anchor,
              use_size = 1, use_shape = 1)
  make_init <- function(chain_id) {
    set.seed(chain_id)
    Z <- Z0 + matrix(rnorm(N * P, sd = 0.05), N, P)
    Z <- sweep(Z, 2, colMeans(Z))
    Z <- Z * sqrt(N / sum(Z^2))
    yr <- matrix(0, N - 1, P)
    v <- numeric(P)
    for (p in seq_len(P)) {
      others <- setdiff(seq_len(N), anchor[p])
      yr[, p] <- Z[others, p]
      v[p] <- log(Z[anchor[p], p] - mean(Z[others, p]))
    }
    dl <- sort(sqrt(rowSums((Z[dat$row, , drop = FALSE] - Z[dat$col, , drop = FALSE])^2)))
    sig <- 0.5
    gam <- seq(dl[2], dl[length(dl) - 1], length.out = dat$K - 1) / sig
    list(yr = yr, v = v, c_raw = rep(0, N), et = replicate(N, rep(0, P), simplify = FALSE),
         m_c = -2, s_c = 0.5, s_e = 0.5, sig = sig, gam = gam)
  }
  list(dat = dat, make_init = make_init, lab = rownames(s))
}

# 書籍 5.1 の Harshman et al. (1982) の再尺度化(横山先生の harshman.R を写したもの) ----------------------
harshman <- function(data, usediag = 0) {
  if (usediag == 0) diag(data) <- NA
  check <- 1
  count <- 0
  const <- 1
  data_sumrc <- rowSums(data, na.rm = TRUE) + colSums(data, na.rm = TRUE)
  oridata_sumrc_ave <- mean(data_sumrc)
  while ((check > 1e-11) && (count < 100)) {
    const_tmp1 <- oridata_sumrc_ave / data_sumrc
    rescaleddata <- data * const_tmp1 %o% const_tmp1
    rescaleddata_sumrc <- rowSums(rescaleddata, na.rm = TRUE) + colSums(rescaleddata, na.rm = TRUE)
    rescaleddata_sumrc_ave <- mean(rescaleddata_sumrc)
    const_tmp2 <- oridata_sumrc_ave / rescaleddata_sumrc_ave
    rescaleddata <- rescaleddata * const_tmp2
    rescaleddata_sumrc <- rowSums(rescaleddata, na.rm = TRUE) + colSums(rescaleddata, na.rm = TRUE)
    check <- sum((rescaleddata_sumrc - data_sumrc)^2)
    count <- count + 1
    data <- rescaleddata
    data_sumrc <- rescaleddata_sumrc
    const <- const * const_tmp1 * sqrt(const_tmp2)
  }
  list(rescaleddata = data, constant = const)
}

# 4つの下位モデルを推定して比べる ----------------------------------------------------------
specs <- tribble(
  ~model,         ~use_size, ~use_shape,
  "楕円(II)",          1,          1,
  "円(III)",           1,          0,
  "等面積の楕円",      0,          1,
  "対称",              0,          0
)
mod <- cmdstanr::cmdstan_model("develop/ellipse_rank.stan", compile_model_methods = TRUE,
                               force_recompile = TRUE)

fit_one <- function(pr, model, use_size, use_shape) {
  od <- tempfile("fit_")
  dir.create(od)
  f <- mod$sample(data = modifyList(pr$dat, list(use_size = use_size, use_shape = use_shape)),
                  init = pr$make_init, chains = 4, parallel_chains = 4, iter_warmup = 1000,
                  iter_sampling = 2500, seed = 1, refresh = 0, show_messages = FALSE,
                  output_dir = od)
  f$init_model_methods()
  U <- unclass(as.matrix(f$unconstrain_draws(format = "draws_matrix")))
  colnames(U) <- paste0("u", seq_len(ncol(U)))
  lb <- setNames(rep(-Inf, ncol(U)), colnames(U))
  ub <- setNames(rep(Inf, ncol(U)), colnames(U))
  lp <- function(x, data) f$log_prob(unconstrained_variables = x, jacobian = TRUE)
  b <- bridge_sampler(samples = U, log_posterior = lp, data = NULL, lb = lb, ub = ub,
                      method = "warp3", repetitions = 5, silent = TRUE)
  N <- pr$dat$N
  phi <- f$draws("phi", format = "matrix")
  ar <- apply(phi, 1, function(v) { m <- matrix(v, N, pr$dat$P); max(apply(m, 1, max) / apply(m, 1, min)) })
  mu <- f$draws("mu", format = "matrix")
  tibble(model = model, logml = mean(b$logml), logml_sd = sd(b$logml),
         divergent = sum(f$diagnostic_summary(quiet = TRUE)$num_divergent),
         s_c = mean(f$draws("s_c")), s_e = mean(f$draws("s_e")), sig = mean(f$draws("sig")),
         max_aspect_med = median(ar), n_mu_neg = sum(colMeans(mu < 0) > 0.5))
}

run_dataset <- function(s, name, type = "similarity") {
  pr <- prep_data(s, type)
  res <- pmap_dfr(specs, function(model, use_size, use_shape) fit_one(pr, model, use_size, use_shape)) %>%
    mutate(data = name, diff_from_best = logml - max(logml), .before = 1)
  res
}

# データ(Rscript で実行したときだけ。source() したときは関数だけを読み込む) ---------------------
if (sys.nframe() == 0) {
args <- commandArgs(trailingOnly = TRUE)
which_data <- if (length(args) > 0) args[1] else "nations"

if (which_data == "nations") {
  # 書籍表1.14: 10か国間の類似度(2ソース)。ソースごとに一相二元として当てる
  x <- read.csv("tests/testthat/fixtures/okada_imaizumi_2024_table1_14.csv", check.names = FALSE)
  out <- map_dfr(split(x, x$source), function(b) {
    m <- as.matrix(b[, -(1:2)])
    rownames(m) <- b$from
    run_dataset(m, paste0("nations_source", b$source[1]))
  })
}
if (which_data == "carswitching") {
  # Harshman et al. (1982) の16車種セグメント間の乗り換え度数(asymmetry パッケージ)。
  # 書籍 5.1 の再尺度化(行 j と列 j に c_j を掛け，[行和 + 列和] をそろえる。対角＝同じセグメントの
  # 買い替えは除く)で市場規模の違いを消してから，類似度として当てる。s_jk / s_kj の比は変わらない
  data(carswitching, package = "asymmetry")
  m <- harshman(as.matrix(carswitching))$rescaleddata
  out <- run_dataset(m, "carswitching_rescaled")
}

if (which_data == "morse") {
  # Rothkopf (1957) のモールス信号の混同(36記号・非類似度・smacof パッケージの morse2)
  data(morse2, package = "smacof")
  out <- run_dataset(as.matrix(morse2), "morse", type = "dissimilarity")
}

saveRDS(out, sprintf("develop/20261121_ellipse_rank_datasets_%s.rds", which_data))
}
