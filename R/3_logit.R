#' Logit transform
#' @description Take the logit transform of a number p
#'
#' @param p Number to logit-transform
#'
#' @example
#' logit(p = 0.01)
#'
#' @export
logit <- function(p){
  r <- log(p/(1 - p))
  return(r)
}
