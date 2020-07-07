#' Get the bycatch mortality rate F that results in a given depletion level
#'
#' This function solves for the bycatch rate \eqn{E} that gives depletion level \code{InitDepl.w}
#'
#' @param f.start An initial guess for the bycatch rate E. The default value is E = 0.5
#' @param S0.w Calf or pup survival, a numeric value between 0 and 1. (Note: the 'w' suffix indicates that \eqn{z} is in the wrapper fn, and is used inside the function by optim())
#' @param S1plus.w  Survival for animals age 1 year and older (a value between 0 and 1)
#' @param K1plus.w The pre-exploitation population size of individuals aged 1 and older. If this value is unavailable, it can be approximated by using the initial depletion and the estimate of current abundance.
#' @param AgeMat.w Age at maturity in years (assumed to be age at first parturition - 1)
#' @param nages.w "Maximum" age, treated as the plus group age. The plus group age can be set equal to the age at maturity +2 years without losing accuracy.
#' @param InitDepl.w The depletion level to solve for (this is equivalent to 'starting depletion')
#' @param z.w The degree of compensation. The default value is \code{z = 2.39}.
#' @param lambdaMax.w The maximum intrinsic growth rate
#' @param N0.w Numbers per recruit in terms of mature individuals individuals
#' @param P0.w  Numbers per recruit in terms of aged 1+ individuals individuals
#'
#' @return A single value, the bycatch rate that would lead to the depletion level \code{InitDepl.w}
#'
#' @examples
#' # Get number of individuals per recruit in terms of mature individuals (\eqn{N0.w})
#' NPROut <- NPR(S0 = .5,S1plus = .944,nages = 25, AgeMat = 18, f=0)
#' N0 <- NPROut$npr # mature numbers per recruit
#' # Get number of individuals per recruit in terms of individuals aged 1+ (\eqn{P0.w})
#' P0 <- NPROut$P1r # 1+ nums per recruit
#'
#' getF(f.start = 0.5,S0.w = 0.5,S1plus.w = 0.944,nages.w = 25,K1plus.w = 9000,AgeMat.w = 18,InitDepl.w = 0.9,z.w = 2.39,lambdaMax.w = 1.04,N0.w = N0,P0.w = P0)
#'
#' @export
getF <- function(f.start = NA, S0.w = NA, S1plus.w = NA, nages.w = NA, K1plus.w = NA, AgeMat.w = NA, InitDepl.w = NA, z.w = NA, lambdaMax.w = NA, N0.w = NA, P0.w = NA){
  # fecundity at unfished equilibrium
  Fec0 = 1.0/N0.w # Fec0=0.012
  FecMax <- getFecmax2(S1plus = S1plus.w,S0 = S0.w,AgeMat = AgeMat.w,lambdaMax = lambdaMax.w)
  A <- (FecMax - Fec0) / Fec0 # Pella-Tomlinson resilience parameter, from Punt 1999 (Annex R) - Note: There are two ways to get A, this is the simpler one

  # set limit for uniroot-- don't search too far in the positive direction
  search.limit <- (1-(1/lambdaMax.w)) * 2.5

  logit.start <- logit(f.start)
  to.minimize <- function(lf = logit.start,
                          S0 = S0.w, S1plus = S1plus.w,
                          nages = nages.w,K1plus = K1plus.w,
                          AgeMat = AgeMat.w, InitDepl = InitDepl.w,
                          z = z.w, lambdaMax = lambdaMax.w,
                          P0 = P0.w, #1+ adults per recruit @E=0
                          N0 = N0.w #mature adults per recruit
  ){
    f <- inv.logit(lf) # f ~= 0.01

    # Get numbers per recruit at f=f
    NE <- NPR(S0 = S0,S1plus = S1plus,nages = nages, AgeMat = AgeMat, f=f)$npr #change this later?
    PE <- NPR(S0 = S0,S1plus = S1plus,nages = nages, AgeMat = AgeMat, f=f)$P1r #new

    R0 = 1 # Can do this because this fn is "scale irrelevant" (AEP)

    # Full solution for R.F is in methods
    if ((1-Fec0*NE)/(Fec0*NE*A) > 1){
      R.F <- 0
    }else{
      R.F <- ( 1 - (1-Fec0*NE)/(Fec0*NE*A))^(1/z) * (R0*P0/PE)}

    # Solve for the R.F that makes Pt1 = Pt2
    Pt1 <- R.F*PE/(R0*P0)
    Pt2 <- InitDepl

    diff <- Pt1-Pt2
    return(diff)
  }

  zero.cross <- uniroot(f = to.minimize,
                        interval = logit(c(0.00001,search.limit)),
                        tol=1e-7)
  f<- inv.logit(zero.cross$root)

  # Check to make sure value of to-minimize is zero
  to.min.val <- to.minimize(lf = logit(f))
  #cat("   val of to.minimize (should be 0)",to.min.val)

  return(f)
}
