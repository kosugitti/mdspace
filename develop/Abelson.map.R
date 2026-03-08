Abelson.map <-
function(dat,locations){
  z <- double()
  X <- dat[,1]
  Y <- dat[,2]
  P <- dat[,3]
  un <- matrix(1,nrow(locations),1,)
  Xs <- un %*% X
  Ys <- un %*% Y
  dm <- ((locations[,1]-Xs)^2+(locations[,2]-Ys)^2)+1
  V <- t(P %*% (1/t(dm)))
  xx <- sort(unique(locations[,1]))
  nx <- length(xx)
  yy <- sort(unique(locations[,2]))
  ny <- length(yy)
  values <- matrix(V,ncol=ny)
  ret <- structure(list(x=xx,y=yy,valence=values))
  return(ret)
}
