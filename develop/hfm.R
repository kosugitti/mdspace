hfm <-
function(data)
{
       if(!is.matrix(data))            #as.matrix
               data <- as.matrix(data)
       if(ncol(data)!=nrow(data))      #check the data whether square matrix or not
               stop("data is not a square matrix")

       #Hermitian form
       Her <- hermitian.form(data)

       #Eigen decomposition
       ev   <- eigen(Her,symmetric=TRUE)
       eval <- ev$values
       evec <- ev$vectors
       rownames(evec) <- colnames(data)
       #Sort by absolute
       od <- abs(eval)
       eval <- eval[order(od,decreasing=TRUE)]
       evec <- evec[,order(od,decreasing=TRUE)]
       #Cal the GOF
       gof   <-eval*eval / (t(eval)%*%eval) 

       return(list(Raw= data, Hermitian =Her ,Eigen=eval,GOF=gof,Vecs=evec))

}
