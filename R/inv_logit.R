# Get inverse logit of a value
#' @description The inverse logit of x
#' @param x A number
#' @return the inverse logit of x. Returns NA if is.na(x)
#' @examples
#' inv.logit(0)
#' inv.logit(1)
inv.logit <- function(x){
  rev <- exp(x)/(1+exp(x))
  return(rev)
}
