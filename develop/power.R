power <-
function(A,max.itr=100,print=FALSE,plot=FALSE,eps=1e-20){

  if(!is.matrix(A))            #as.matrix
    A <- as.matrix(A)
  if(ncol(A)!=nrow(A))      #check the data whether square matrix or not
    stop("data is not a square matrix")

  u <- list(rep(1,length=nrow(A)),max.itr)
  v <- list(rep(0,length=nrow(A)),max.itr)
  e.val <- 0
  for(i in 1:max.itr){
    v[[i]]<- A%*%u[[i]]
    if(abs(max(v[[i]])-e.val)<eps){break;}
    e.val <- max(v[[i]])
    u[[i+1]]<- v[[i]]/e.val
    e.vec <- u[[i+1]]
    if(print){cat("itr",i,":",sprintf("%10.10f",e.val),"\n")}
  }
  if(plot){
    # chaning eigenvector matrix
    change <- matrix(unlist(v),nrow=nrow(A))
    # plotting eigenvalue of each itterations
    plot(apply(change,2,max),type="b")
  }
  
  res <- list(E.val = e.val,E.vec=e.vec,n.itr=i)
  return(res)
  
}
