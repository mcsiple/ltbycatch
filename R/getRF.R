#' Get recruitment at exploitation rate E
#'
#' @description calculates recruitment at exploitation rate E
#' @param E initial exploitation rate (num)
#' @param S0 Calf or pup survival (num)
#' @param S1plus Adult survival  (num)
#' @param nages Number of age classes, incl plus group age (num)
#' @param K1plus Adult carrying capacity (num)
#' @param AgeMat Age at maturity (num)
#' @param z degree of compensation (num)
#' @param A Pella-Tomlinson resilience parameter (Punt 1999; Annex R). A = (FecMax - Fec0) / Fec0 (num)
#' @param P0 Unfished 1+ numbers per recruit, \tildeP(0). Must double-check with AEP.
#' @param N0 Unfished mature numbers per recruit, \tildeN(0). Must double-check with AEP.
#' @return recruitment given exploitation rate E - to multiply by Ninit to get initial nums at age (num)
#'
#' @export
getRF <- function(E, S0, S1plus, nages,K1plus, AgeMat, z, A, P0, N0){

  NE <- NPR(S0 = S0,S1plus = S1plus,nages = nages, AgeMat = AgeMat, f=E)$npr
  PE <- NPR(S0 = S0,S1plus = S1plus,nages = nages, AgeMat = AgeMat, f=E)$P1r

  # More general version ('actual' R0 = K1plus/P1r)
  R0 <- 1

  # Fecundity at unfished equilibrium
  Fec0 = 1.0/N0

  if ((1-Fec0*NE)/(Fec0*NE*A) > 1)
    R.F <- 0
  else
    R.F <- ( 1 - (1-Fec0*NE)/(Fec0*NE*A))^(1/z) * ((R0*P0)/PE) # this used to be NE!

  return(R.F)
}
