inv.logit <- function(x){
  #' @description takes the inverse logit of x
  #' @param x number to take inv logit of (num)
  #' @return inverse logit of x, returns NA if is.na(x)
  rev <- exp(x)/(1+exp(x))
  return(rev)
}
