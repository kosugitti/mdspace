#' Dual Sacling with EigenDecomposition
#'
#' @param dat A matrix or data frame
#' @export


DualScaling <- function(dat) {
  # 与えられたデータの隔離
  nc <- ncol(dat)
  nr <- nrow(dat)
  cname <- colnames(dat)
  rname <- rownames(dat)
  # 抽出する次元数
  dm <- min(nc, nr) - 1
  # 出力するベクトルの準備
  normed.col <- matrix(ncol = nc, nrow = dm)
  projected.col <- matrix(ncol = nc, nrow = dm)
  normed.row <- matrix(ncol = nr, nrow = dm)
  projected.row <- matrix(ncol = nr, nrow = dm)
  singular <- c()

  tmp <- dat
  # 標準形からの固有値分解
  D <- diag(apply(tmp, 2, sum))
  Dn <- diag(apply(tmp, 1, sum))
  Dhlf <- sqrt(D)
  Dm_hlf <- diag(diag(1 / Dhlf))
  C0 <- Dm_hlf %*% t(tmp) %*% solve(Dn) %*% tmp %*% Dm_hlf
  C1 <- C0 - (Dhlf %*% matrix(1, nrow = nrow(tmp), ncol = ncol(tmp)) %*% Dhlf) / sum(tmp)
  EigenSystem <- eigen(C1)
  # 有効次元までの固有値抽出
  eigenvalue <- EigenSystem$values[1:dm]
  singular <- sqrt(eigenvalue)
  # 固有ベクトルの抽出
  vec <- EigenSystem$vectors[, 1:dm] * sqrt(sum(tmp))
  # 列スコア
  nu <- Dm_hlf %*% vec
  # 投影された行スコア
  pv <- solve(Dn) %*% tmp %*% nu
  # 行スコア
  nv <- pv / rep(1, nr) %*% t(singular)
  # 投影された列スコア
  pu <- solve(D) %*% t(tmp) %*% nv

  # Correct dims
  delta_k <- singular^2 / sum(singular^2) * 100
  cum_delta_k <- c()
  cum_delta_k[1] <- delta_k[1]
  for (i in 2:dm) {
    cum_delta_k[i] <- cum_delta_k[i - 1] + delta_k[i]
  }

  # var names
  colnames(normed.col) <- cname
  colnames(projected.col) <- cname
  colnames(normed.row) <- rname
  colnames(projected.row) <- rname
  dimName <- paste0("Dim", 1:dm)
  rownames(normed.row) <- dimName
  rownames(normed.col) <- dimName
  rownames(projected.row) <- dimName
  rownames(projected.col) <- dimName

  # Style the output
  res <- t(cbind(singular, singular^2, delta_k, cum_delta_k))
  colnames(res) <- dimName
  rownames(res) <- c("singular value", "eigen value", "contribution", "cumulative contribution")

  normed.col <- t(normed.col)
  normed.row <- t(normed.row)
  projected.col <- t(projected.col)
  projected.row <- t(projected.row)



  return(list(
    result = res,
    ndims = dm,
    singularValue = singular,
    eigenValue = singular^2,
    delta = delta_k,
    cumulativeDelta = cum_delta_k,
    NormedCol = nu,
    NormedRow = nv,
    ProjectedCol = pu,
    ProjectedRow = pv
  ))
}
