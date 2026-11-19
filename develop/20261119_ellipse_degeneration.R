# ellipse_mds() の試作で見つかった細長い楕円への流れ(2026-11-19) ------------------------------
# (1) kappa(縦横比の上限)を変えたときの S と半軸の大きさ (2) 本の表4.4 の解から最小化し直したときの行き先
devtools::load_all(quiet = TRUE)
occ <- as.matrix(read.csv("tests/testthat/fixtures/okada_imaizumi_2024_table1_19.csv", row.names = 1, check.names = FALSE))
for (k in c(3, 5, 10, 20, Inf)) {
  tm <- system.time(f <- suppressWarnings(ellipse_mds(occ, kappa = k, nrandom = 30, check_orientation = FALSE)))
  sa <- f$semiaxis
  cat(sprintf("kappa=%s S=%.4f at_bound=%s max_ratio=%.1f max_semiaxis=%.2f n<=.3245=%d time=%.1f\n", k, f$stress, f$at_bound,
              max(apply(sa, 1, max) / apply(sa, 1, min)), max(sa), sum(f$starts <= 0.3245), tm["elapsed"]))
}
occ <- as.matrix(read.csv("tests/testthat/fixtures/okada_imaizumi_2024_table1_19.csv", row.names = 1, check.names = FALSE))
tab <- matrix(c(1.191,0.656,1.665,0.624, 1.025,0.047,1.290,0.394, 0.988,-1.020,1.165,1.154,
 -0.188,0.395,1.143,0.599, -0.311,-0.342,0.149,0.175, 0.215,-0.689,0.262,0.117,
 -0.453,0.956,1.175,0.821, -1.241,0.355,1.072,1.303, -0.703,-0.697,0.679,0.430,
 -0.523,0.339,1.091,0.628), 10, byrow = TRUE)
use <- row(occ) != col(occ)
cells <- list(s = occ[use], j = row(occ)[use], k = col(occ)[use], n = 10, decreasing = TRUE,
              model = "II", alpha = -1, beta = 1, kappa = Inf)
f <- mdspace:::.ell_optim(tab[, 1:2], tab[, 3:4], cells, 5000, 1e-12)
nr <- mdspace:::.ell_normalize(f$conf, f$axis, "II")
cat("from book: S", f$stress, "\n"); print(round(cbind(nr$conf, nr$axis), 3))
