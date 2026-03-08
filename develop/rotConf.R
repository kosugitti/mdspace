rotConf <-
function(map,dig){
  rad <- -dig*pi/180
  rotmat <- matrix(c(cos(rad),-sin(rad),sin(rad),cos(rad)),byrow=T,2,2)
  R_map <- map %*% rotmat
  return(R_map)
}
