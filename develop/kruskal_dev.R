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
mat.dist <- as.dist(mat)


### initialize
nobs <- nrow(mat)
max_dim <- 3
min_dim <- 3
max_iter <- 100
eps <- 1e-3
eps_diff <- 1e-7
n_dim <- max_dim

### Initial configuration
conf <- if (n_dim == max_dim) {
    cmdscale(mat.dist, k = n_dim)
} else {
    tmp <- prcomp(mat.dist)
    tmp$x[, 1:n_dim]
}

alpha <- 0.2
stress <- numeric(max_iter)
FLG <- TRUE
iter <- 1

# Precompute values for tmp.df creation
i_vals <- rep(1:(nobs - 1), sapply(1:(nobs - 1), function(x) nobs - x))
j_vals <- unlist(sapply(1:(nobs - 1), function(x) (x + 1):nobs))

# Main loop
while (FLG) {
    # Calculate the distance matrix
    n_dist <- as.matrix(dist(conf))

    # Efficient creation of tmp.df
    oij_vals <- mapply(function(i, j) mat[j, i], i_vals, j_vals)
    dij_vals <- mapply(function(i, j) n_dist[j, i], i_vals, j_vals)
    tmp.df <- data.frame(i = i_vals, j = j_vals, oij = oij_vals, dij = dij_vals)

    # Disparity calculation
    tmp.df <- tmp.df %>% arrange(oij, dij)
    avg_dij <- tapply(tmp.df$dij, tmp.df$oij, mean)

    # Adjust avg_dij
    for (i in 2:length(avg_dij)) {
        if (avg_dij[i] < avg_dij[i - 1]) {
            new_avg <- mean(c(avg_dij[i], avg_dij[i - 1]))
            avg_dij[i] <- new_avg
            avg_dij[i - 1] <- new_avg
            if (i > 2) {
                i <- i - 1
            } else {
                i <- i + 1
            }
        }
    }

    # Map adjusted values back to tmp.df
    tmp.df$disp <- avg_dij[tmp.df$oij]

    # Stress calculation
    num <- sum((tmp.df$dij - tmp.df$disp)^2)
    den <- sum((tmp.df$dij - mean(tmp.df$dij))^2)
    S <- stress[iter] <- sqrt(num / den)
    g <- matrix(0, nrow = nrow(conf), ncol = ncol(conf))

    # Update configurations
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
                        A <- (dij.tmp - disp.tmp) / num
                        B <- dij.tmp / den
                        C <- (conf[i, l] - conf[j, l]) / dij.tmp
                        tmp <- tmp + (Kron * (A - B) * C)
                    }
                }
            }
            g[k, l] <- S * tmp * -1
        }
    }

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


    # Judge termination
    if (S < eps) {
        print("Stress < eps!")
        FLG <- FALSE
    } else if (iter > 1 && abs(stress[iter] - stress[iter - 1]) < eps_diff) {
        print("Stress change < eps!")
        FLG <- FALSE
    } else if (iter == max_iter) {
        print("reach iter max")
        FLG <- FALSE
    }

    iter <- iter + 1
}

