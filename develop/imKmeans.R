#' improved k-means function
#' 
#' @param dat data as matrix
#' @param cl number of clusters
#' @return cluster cluster that assigned for each obs.
#' @return aic AIC
#' @return bic BIC
#' 
#' @export
#' @examples 
#' imKmeans(iris[1:4],3)

imKmeans <- function(dat,cl){
  # step 1
  dat.km <- kmeans(dat,cl,algorithm = "MacQueen")
  # step 2
  conv <- 3
  clust <- dat.km$cluster
  
  while(conv>0){
    subdat <- list()
    mu <- list()
    Sig <- list()
    kyori <- list()
    lk <- list()
    for(i in 1:cl){
      subdat[[i]] <- dat[clust==i,]
      mu[[i]] <- apply(subdat[[i]],2,mean)
      Sig[[i]] <- cov(subdat[[i]])
      #distance and liklihood
      dat.m <- as.matrix(dat)
      kyori[[i]] <- mahalanobis(dat,mu[[i]],Sig[[i]])
      lk[[i]] <- apply(dat.m,1,function(j){ #likelihood
        (2 * pi)^(-nrow(dat)/2) * det(Sig[[i]])^(-0.5) * exp(-0.5*(t(j-mu[[i]]) %*% solve(Sig[[i]],LINPACK=TRUE) %*% (j-mu[[i]])))
      })
    }
    lk.df <- as.data.frame(matrix(unlist(lk),ncol=cl))
    # step3
    neo.clust <- apply(lk.df,1,which.max)
    conv <- conv + ifelse(sum(abs(neo.clust-clust))!=0,1,0)
    if(sum(abs(neo.clust-clust))==0){
      conv <- conv - 1
      if(conv==0){break}
    }else{
      clust<-neo.clust
    }
  }
  
  # step.4
  clust.ratio <- table(clust)/nrow(dat)
  
  # step.5 probability dens.
  Ld <- list()
  for(i in 1:cl){
    Ld[[i]] <- apply(dat,1,function(j){ #prob.density
      mvtnorm::dmvnorm(j,mu[[i]],Sig[[i]])*clust.ratio[i]
    })
  }
  L <- apply(as.data.frame(matrix(unlist(Ld),ncol=cl)),1,sum) # lk
  l <- (ncol(dat)*(ncol(dat)+3))*cl*0.5
  (aic <- -2* sum(log(L)) + 2*l)
  (bic <- -2* sum(log(L)) + l*(log(l)))
  
  return(list(cluster=clust,aic=aic,bic=bic))
  
}