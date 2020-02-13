#' Logit transform
#' @description Take the logit transform of a number p
#' @param p Number to take the logit of
#'
#' @export
logit <- function(p){
  r <- log(p/(1 - p))
  return(r)
}
