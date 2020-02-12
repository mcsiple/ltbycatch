#' Get the value of F that results in a given depletion level
#'
#' @description solves for f that gives depletion level InitDepl.w
#' @param f.start initial guess for F (num)
#' @param S0.w Calf or pup survival. 'w' suffix indicates that z is in the wrapper fn, and is used inside the fn by optim()  (num)
#' @param S1plus.w Adult survival (num)
#' @param K1plus.w Adult carrying capacity (num)
#' @param AgeMat.w Age at maturity (generally AFP+1) (num)
#' @param InitDepl.w Depletion level to solve for (num)
#' @param z.w degree of compensation (num)
#' @param lambdaMax.w max intrinsic growth rate (num)
#' @param unpr.w unfished nums per recruit (num)

getF <- function(f.start, S0.w, S1plus.w, nages.w,K1plus.w, AgeMat.w, InitDepl.w, z.w,lambdaMax.w, unpr.w){
  # Calcs now outside this function: unfished nums per recruit

  # fecundity at unfished equilibrium
  Fec0 = 1.0/unpr.w # Fec0=0.012
  FecMax <- getFecmax2(S1plus = S1plus.w,S0 = S0.w,AgeMat = AgeMat.w,lambdaMax = lambdaMax.w)
  A <- (FecMax - Fec0) / Fec0 # Pella-Tomlinson resilience parameter, from Punt 1999 (Annex R) - Note: There are two ways to get A, this is the simpler one

  # set limit for uniroot-- don't search too far in the positive direction
  search.limit <- (1-(1/lambdaMax.w)) * 2.5

  logit.start <- logit(f.start)
  to.minimize <- function(lf = logit.start,
                          S0 = S0.w, S1plus = S1plus.w,
                          nages = nages.w,K1plus = K1plus.w,
                          AgeMat = AgeMat.w, InitDepl = InitDepl.w,
                          z = z.w, lambdaMax = lambdaMax.w, unpr = unpr.w){
    f <- inv.logit(lf) # f ~= 0.01

    # Get numbers per recruit at f=f
    nprf <- NPR(S0 = S0,S1plus = S1plus,nages = nages, AgeMat = AgeMat, f=f)$npr
    R0 <- 1 # Can do this because this fn is "scale irrelevant" (AEP)

    # Full solution for R.F is in methods

    if ((1-Fec0*nprf)/(Fec0*nprf*A) > 1){
      R.F <- 0
    }else{
      R.F <- ( 1 - (1-Fec0*nprf)/(Fec0*nprf*A))^(1/z) * ((R0*unpr)/nprf)}

    # Solve for the R.F that makes the following statement true
    Pt1 <- R.F*nprf/(R0*unpr)
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
