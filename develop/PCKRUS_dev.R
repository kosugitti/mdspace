rm(list = ls())
pacman::p_load(tidyverse, MASS)
dat <- c(3, 5, 1, 6, 5, 6, 7, 4, 1, 2, 5, 2, 4, 6, 6, 1, 3, 2, 7, 4, 3)
name <- c("Greece", "Hawai", "WestCoast", "HongKong", "London", "EastCoast", "Australia")
nc <- 7

dat <- c(7.53319,8.12766,5.29787,5.56383,2.81915,6.14894,2.17021,6.64894,3.87234,6.12766)
nc <- 5


mat <- matrix(0, ncol = nc, nrow = nc)

# dat <- c(3,2,5,1,4,6)
# mat <- matrix(0,4,4)
# nc <- 4
#
pos <- 0
for(i in 1:(nc-1)){
    for(j in (i+1):nc){
        pos <- pos + 1
        mat[j,i] <- dat[pos]
    }
}
# colnames(mat) <- rownames(mat) <- name
# mat <- mat + t(mat)
mat.dist <- as.dist(mat)

# result.iso <- MASS::isoMDS(d=mat.dist,k=2)

# mat <- as.matrix(dist(as.matrix(swiss[, -1])))
# mat.dist <- as.dist(mat)
#
# mat[lower.tri(mat)] <- dat
# mat <- pmax(mat,t(mat))

### initialize
nobs <- nrow(mat)
max_dim <- 3
min_dim <- 3
max_iter <- 100
eps <- 1e-3
eps_diff <- 1e-7
n_dim <- max_dim

### Loop
# iniital configuration
if (n_dim == max_dim) {
  conf <- cmdscale(mat.dist, k = n_dim)
} else {
  tmp <- prcomp(mat.dist)
  conf <- tmp$x[, 1:n_dim]
}

# conf <- matrix(c(3,2,1,10,2,7,3,4),ncol=2)

alpha <- 0.2
stress <- rep(NA, max_iter)
FLG <- TRUE
iter <- 1

while (FLG) {
  # calc the distance
  n_dist <- as.matrix(dist(conf))
  nr <- nobs * (nobs-1) /2
  tmp.df <- data.frame(origin = 1:nr,i=NA,j=NA,oij=NA,dij=NA)
  pos <- 0
  for(i in 1:(nobs-1)){
      for(j in (i+1):nobs){
          pos <- pos + 1
          tmp.df[pos,]$i <- i
          tmp.df[pos,]$j <- j
          tmp.df[pos,]$oij <- mat[j,i]
          tmp.df[pos,]$dij <- n_dist[j,i]
      }
  }

  # disparity
  # oij <- 1:10
  # dij <- c(2,4,2,5,7,5,2,9,10,6)
  # oij と dij を組み合わせたデータフレームを作成
  tmp.df <- tmp.df %>% arrange(oij, dij)
  # ソートされたデータで avg_dij を計算
  avg_dij <- tapply(tmp.df$dij, tmp.df$oij, mean)

  # oij の順序に従って avg_dij を調整するための反復処理
  i <- 2
  while (i <= length(avg_dij)) {
    if (avg_dij[i] < avg_dij[i - 1]) {
      new_avg <- mean(c(avg_dij[i], avg_dij[i - 1]))
      avg_dij[i] <- new_avg
      avg_dij[i - 1] <- new_avg
      if (i > 2) {
        i <- i - 1 # 一つ前に戻る
      } else {
        i <- i + 1
      }
    } else {
      i <- i + 1
    }
  }

  # 調整された平均値を元のデータフレームにマッピング
  tmp.df$disp <- avg_dij[tmp.df$oij]

  # print(tmp.df)
  # plot(tmp.df$disp,tmp.df$oij)

  # calc the stress
  num <- (tmp.df$dij - tmp.df$disp)^2 |> sum()
  # stress2
  den <- (tmp.df$dij - mean(tmp.df$dij))^2 |> sum()
  # stress1
  # den <- tmp.df$dij^2 |> sum()
  S <- stress[iter] <- sqrt(num / den)

  # improve the configurations

  S_star <- num
  # T_star <- tmp.df$dij^2 |> sum()
  T_star <- den

  g <- matrix(0, nrow = nrow(conf), ncol = ncol(conf))

  # webの方法
  # for(i in 1:nobs){
  #     for(l in 1:n_dim){
  #         tmp <- 0
  #         for(j in 1:nobs){
  #             if(j!=i){
  #                 if(i>j){
  #                     ii <- j
  #                     jj <- i
  #                 }else{
  #                     ii <- i
  #                     jj <- j
  #                 }
  #                 dij.tmp <- tmp.df[tmp.df$i==ii & tmp.df$j==jj,]$dij
  #                 disp.tmp<- tmp.df[tmp.df$i==ii & tmp.df$j==jj,]$disp
  #                 tmp <- tmp + (1-disp.tmp/dij.tmp) *(conf[j,l]-conf[i,l])
  #             }
  #         }
  #         g[i,l] <- tmp
  #     }
  # }

  # Kruskalの論文・ただし勾配は逆転
  for (k in 1:nobs) {
    for (l in 1:n_dim) {
      tmp <- 0
      for (i in 1:nobs) {
        for (j in 1:nobs) {
          if (i != j) {
            Kronecker_ki <- ifelse(k == i, 1, 0)
            Kronecker_kj <- ifelse(k == j, 1, 0)
            Kron <- Kronecker_ki - Kronecker_kj
            signum <- sign(conf[i, l] - conf[j, l])
            if(i>j){
                    ii <- j
                    jj <- i
            }else{
                    ii <- i
                    jj <- j
            }
            dij.tmp <- tmp.df[tmp.df$i==ii & tmp.df$j==jj,]$dij
            disp.tmp<- tmp.df[tmp.df$i==ii & tmp.df$j==jj,]$disp
            A <- (dij.tmp - disp.tmp) / S_star
            B <- dij.tmp / T_star
            C <- (conf[i, l] - conf[j, l]) / dij.tmp
            tmp <- tmp + (Kron * (A - B) * C)
          }
        }
      }
      g[k, l] <- S * tmp * -1
    }
  }

  #　高橋の方法
  # for(i in 1:nobs){
  #     for(l in 1:n_dim){
  #         tmp <- 0
  #         for(k in 1:nobs){
  #             if(k!=i){
  #                 if(i>k){
  #                     ii <- k
  #                     jj <- i
  #                 }else{
  #                   ii <- i
  #                   jj <- k
  #                 }
  #               dij.tmp <- tmp.df[tmp.df$i==ii & tmp.df$j==jj,]$dij
  #               disp.tmp<- tmp.df[tmp.df$i==ii & tmp.df$j==jj,]$disp
  #               A <- (dij.tmp - disp.tmp)/S_star
  #               B <- (dij.tmp - mean(tmp.df$dij))/T_star
  #               signum <- sign(conf[i,l]-conf[k,l])
  #               tmp <- tmp + signum * ((1/dij.tmp) * (A-B) * abs(conf[i,l]-conf[k,l]))
  #             }
  #         }
  #         g[i,l] <- 2 * S * tmp * -1
  #     }
  # }

  mag <- sqrt(sum(g^2)) / sqrt(sum(conf^2))


  # angle * relax * goodluck
  if (iter == 1) {
    angle <- 1
  } else {
    angle <- sum(g * old_g) / (sqrt(sum(g^2)) * sqrt(sum(old_g^2)))
  }
  old_g <- g

  if (iter < 5) {
    S_5 <- S
    relax <- min(1, (S / S_5))
  } else {
    relax <- min(1, (S / stress[iter - 5]))
  }

  if (iter == 1) {
    goodluck <- 1
  } else {
    goodluck <- min(1, (S / stress[iter - 1]))
  }

  alpha <- alpha * angle * relax * goodluck

  # update
  conf <- conf + (alpha/mag) * g


  # judge -----------------------------------------------------------

  if (S < eps) {
      print("Stress < eps!")
    FLG <- FALSE
  }
  if (iter > 1) {
    if (abs(stress[iter] - stress[iter - 1]) < eps_diff) {
        print("Stress change < eps!")
      FLG <- FALSE
    }
  }
  if (iter == max_iter) {
      print("reach iter max")
    FLG <- FALSE
  }

  iter <- iter + 1
  print(S)
}


plot(stress)
resultISO <- MASS::isoMDS(as.dist(mat), k = 2)
plot(resultISO$points[,1],resultISO$points[,2])
plot(-conf[,1],-conf[,2])
