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
#' @param unpr Unfished numbers-per-recruit (num)
#' @return recruitment given exploitation rate E - to multiply by Ninit to get initial nums at age (num)
#'
#' @export
getRF <- function(E, S0, S1plus, nages,K1plus, AgeMat, z, A, unpr){ #I have remvoed lambdaMax because it's only used here to calculate A, which is now done outside the fn!

  nprf <- NPR(S0 = S0,S1plus = S1plus,nages = nages, AgeMat = AgeMat, f=E)$npr

  # More general version ('actual' R0 = K1plus/P1r)
  R0 <- 1

  # Fecundity at unfished equilibrium
  # AEP: You could pass FecMax and A into the function because all you are changing is E
  # MCS: I don't use FecMax except to calculate A so the fn now just passes A
  Fec0 = 1.0/unpr

  if ((1-Fec0*nprf)/(Fec0*nprf*A) > 1)
    R.F <- 0
  else
    R.F <- ( 1 - (1-Fec0*nprf)/(Fec0*nprf*A))^(1/z) * ((R0*unpr)/nprf)

  return(R.F)
}
