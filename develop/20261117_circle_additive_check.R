# circle_mds の初期布置: Torgerson の加算定数なし/ありの比較(合理的初期値だけ・nrandom = 0 相当)
suppressMessages(pkgload::load_all(".", quiet = TRUE))
# 合理的初期値だけ(nrandom = 0)で，加算定数なし/ありを比べる
fit_rat <- function(S, p, additive) {
  X <- as_proximity(S, type = "similarity")
  use <- row(X) != col(X) & !is.na(X)
  cells <- list(s = X[use], j = row(X)[use], k = col(X)[use], n = nrow(X), decreasing = TRUE)
  init <- mdspace:::.circle_init(X, p, "similarity", additive = additive)
  mdspace:::.circle_optim(init$conf, init$radius, cells, 1000, 1e-10)$stress
}
occ <- as.matrix(read.csv("tests/testthat/fixtures/okada_imaizumi_2024_table1_19.csv", row.names = 1, check.names = FALSE))
occ_res <- sapply(1:5, function(p) c(p = p, none = fit_rat(occ, p, FALSE), add = fit_rat(occ, p, TRUE)))
conds <- expand.grid(n = c(5, 10, 20), p = c(2, 3, 5, 10), noise = c(0, 0.1)); conds <- conds[conds$p < conds$n, ]
res <- do.call(rbind, lapply(seq_len(nrow(conds)), function(k) {
  n <- conds$n[k]; p <- conds$p[k]; noise <- conds$noise[k]
  out <- t(sapply(1:20, function(i) {
    set.seed(1000 * n + 100 * p + 10 * noise * 10 + i)
    X <- matrix(rnorm(n * p), n, p); D <- as.matrix(dist(X))
    r <- runif(n, 0, 0.3 * median(D[upper.tri(D)])); M <- D - outer(r, r, "-")
    off <- row(M) != col(M); E <- matrix(rnorm(n * n, 0, noise * sd(M[off])), n, n)
    S <- -(M + E); diag(S) <- NA
    use <- row(S) != col(S)
    cells <- list(s = S[use], j = row(S)[use], k = col(S)[use], n = n, decreasing = TRUE)
    st <- mdspace:::.circle_stress(c(X, r), cells, p)
    c(none = fit_rat(S, p, FALSE) > st + 1e-4, add = fit_rat(S, p, TRUE) > st + 1e-4)
  }))
  data.frame(n = n, p = p, noise = noise, local_none = sum(out[, 1]), local_add = sum(out[, 2]))
}))
list(occupation = round(t(occ_res), 4), recovery_rational_only = res, total = colSums(res[, 4:5]))
