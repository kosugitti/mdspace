#' make dominance
#'
#' @param dat A matrix or data frame
#' @export
#'
make.dominance <- function(dat) {
  nc <- ncol(dat)
  nr <- nrow(dat)
  dom <- matrix(ncol = ncol(dat), nrow = nrow(dat))
  for (i in 1:nc) {
    for (j in 1:nr) {
      dom[j, i] <- nc + 1 - 2 * dat[j, i]
    }
  }
  return(dom)
}
