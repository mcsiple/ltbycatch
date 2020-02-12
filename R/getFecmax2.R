# Get max theoretical fecundity -------------------------------------------
# From a ginormous bunch of algebra we did at OMF2 - IMPORTANT: when applying this calculation, use AgeMat and when defining probability of giving birth you use parturition ages (AgeMat+1)

getFecmax2 <- function(S0,lambdaMax,S1plus,AgeMat){
  fmax = (lambdaMax^(AgeMat)-(S1plus*(lambdaMax^(AgeMat-1)))) / (S0*S1plus^(AgeMat-1))
  return(fmax)
}
