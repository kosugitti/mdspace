#' make dummy
#'
#' @param dat A matrix or data frame
#' @export

make.dummy <- function(dat) {
  NC <- NCOL(dat)
  new.nc <- apply(dat, 2, max)
  start.point <- c(rep(-new.nc[1], NC))
  for (i in 1:NC) {
    for (j in 1:i) {
      start.point[i] <- start.point[i] + new.nc[j]
    }
  }
  new.dat <- matrix(0, nrow = NROW(dat), ncol = new.nc * NC)
  for (i in 1:nrow(dat)) {
    for (j in 1:NC) {
      new.dat[i, start.point[j] + dat[i, j]] <- 1
    }
  }
  return(new.dat)
}
