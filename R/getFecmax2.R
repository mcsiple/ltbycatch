#' Get max theoretical fecundity
#'
#' @description Calculate maximum theoretical fecundity Fmax
#' @details From a ginormous bunch of algebra we did at OMF2 - IMPORTANT: when applying this calculation, use AgeMat and when defining probability of giving birth you use parturition ages (AgeMat+1)
#' @param S0 Calf survival (num)
#' @param lambdaMax Maximum population growth rate (num)
#' @param S1plus Survival of age 1+ individuals (num)
#' @param AgeMat Age at maturity (num)
#' @return A list of numbers per recruit (npr), 1+ numbers per recruit (P1r), and numbers at age per recruit (nvec)
#'
#' @example
#' x <- getFecmax2(lambdaMax=1.00,S0=0.944,S1plus=0.99,AgeMat=18)
#' print(x)
#' unpr <- NPR(S0 = 0.944,S1plus = 0.99,AgeMat=18, nages = 10000,f=0)
#' print(1/unpr$npr)
#' @export
getFecmax2 <- function(S0,lambdaMax,S1plus,AgeMat){
  fmax = (lambdaMax^(AgeMat)-(S1plus*(lambdaMax^(AgeMat-1)))) / (S0*S1plus^(AgeMat-1))
  return(fmax)
}
