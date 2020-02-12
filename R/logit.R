logit <- function(p){
  #' @description logit transformation
  r <- log(p/(1 - p))
  return(r)
}
