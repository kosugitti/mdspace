hfm.plot <-
function(data,dim,Xlim=c(-1,1),Ylim=c(-1,1))
{
	plot(data$Vecs[,dim],
 	      main=paste("Dim", dim, "with Eigenvalue" ,round(data$Eigen[dim],3)),
	      xlab="real",    #X-Label
	      ylab="imag",  #Y-Label
	      xlim=Xlim,  #X-range
	      ylim=Ylim,  #Y-range
       	axes=F)
	axis(1, pos = 0, at = -3:3, adj = 0, col = 1) #  Draw X axsis
	axis(2, pos = 0, at = -3:3, adj = 1, las = 2) #  Draw Y axsis
	for( i in seq(along = data$Vecs[,dim]))
		  arrows(0,0,Re(data$Vecs[i,dim]),Im(data$Vecs[i,dim]))
        text(data$Vecs[,dim],rownames(data$Vecs))
}
