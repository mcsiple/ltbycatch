#' Derivative of z
#'
#' @param z
#' @param MNPL
#' @param lh.params
#'
#' @return
#' @export
#'
#' @examples
get_dz <- function(z, MNPL, lh.params){
  # want diff between MNPL and f(z) equal to zero
  lh.params$z = z
  MNPL_calc <- getMNPL(lh.params = lh.params)
  dZ <- MNPL - MNPL_calc
  return(dZ)
}
