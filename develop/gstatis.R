gstatis <-
function(data,size,rep,viewpoint)
{
  # Make data into List form as repeated "rep" times size*size matrix
  # Confirm the argumented values
  if(nrow(data)!=size*rep)
      stop("data size error.Ncol is not match by size*rep")
  
  data.t <- t(data)
  dim(data.t)<-c(size,size,rep)
  gMatrix <- list()
    for(i in 1:rep){
      gMatrix[[i]] <- hermitian.form(as.matrix(data.t[,,i]))
    }
  
  #step 1;PCA of GeneralizedRV 
  GRV <- matrix(ncol=rep,nrow=rep)
  for( i in 1:rep ){
    for( j in 1:rep ){
      GRV[i,j] <- mat.trace(gMatrix[[i]],gMatrix[[j]])/((mat.trace(gMatrix[[i]],gMatrix[[i]])*mat.trace(gMatrix[[j]],gMatrix[[j]]))^0.5)
    }
  }

  GRV.eStructure <- eigen(GRV)
  
  #step 2;intra-structures with point of view
  cH.Matrix <- matrix(0,ncol=size,nrow=size)
  for(i in 1:rep){
    cH.Matrix <- cH.Matrix + (GRV.eStructure$vectors[i,viewpoint] * gMatrix[[i]])
  }
  cH.eStructure <- eigen(cH.Matrix)
  
  #step 3;trajactory analysis
  Opa <- cH.eStructure$vectors
  Opa_ct <- t(Conj(Opa))  #conjugate transpose of O_p^(a)
  Psy <- matrix(diag(cH.eStructure$values),nrow=size)
  Q_a <- Opa %*% solve(Psy) %*% Opa_ct  # orthogonal projector

  tMatrix<-list()
  for(i in 1:rep){
    tMatrix[[i]] <-Q_a %*%t(gMatrix[[i]])
  }
  
  return(list(Hermitian.Form=gMatrix,GRV=GRV,Step1.eigenStructure=GRV.eStructure,           #result of Step1
              constructed.Matrix=cH.Matrix,Step2.eigenStructure=cH.eStructure,              #result of Step2
              tMatrix=tMatrix                                                               #result of Step3
              ))
}
