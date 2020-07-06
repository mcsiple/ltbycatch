#' Generate one marine mammal population trajectory
#'
#' This function generates one trajectory for a marine mammal population, starting at a user-specified depletion level `InitDepl`.
#' @authors
#' @details
#' The population model is a single-sex age-structured model in which the number of calves or pups born each year is density dependent, with the extent of density dependence a function of the number of mature adults \eqn{\tildeN}, the fecundity (pregnancy rate) at pre-exploitation equilibrium \eqn{f_0}, the maximum theoretical fecundity rate fmax, the degree of compensation \eqn{z}, and the abundance of individuals aged 1+ \eqn{N_{y+1}^{1+}} relative to carrying capacity \eqn{K^{1+}}.
#'
#' @param S0 Calf survival rate
#' @param S1plus Survival rate for animals 1 and older
#' @param K1plus The pre-exploitation population size (1+ component of pop)
#' @param AgeMat Age at maturity
#' @param InitDepl Starting depletion level
#' @param ConstantCatch Total bycatch each year, expressed as a vector of length [nyears]
#' @param ConstantF vector (length=nyears) rate of bycatch each year
#' @param z Degree of compensation (2.39 for bowhead whales, set as default for others)
#' @param nyears Number of years to project
#' @param nages "Maximum" age, treated as the plus group age. The plus group age can be set equal to the age at maturity + 2 without losing accuracy.
#' @param lambdaMax max intrinsic growth rate
#'
#' @return a list containing a matrix N of numbers at age (dimensions nyears (rows) x nages (columns)) and one vector TotalPop (a vector of length nyears), containing the number of age 1+ individuals in the population.
#'
#' @examples
#' Dynamics()
#'
#' @export
Dynamics <- function(S0 = NA, S1plus = NA, K1plus = NA, AgeMat = NA, InitDepl = NA, ConstantCatch = NA, ConstantF = NA, z = 2.39, nyears = NA, nages = NA, lambdaMax = NA){
  # Checks
  if(length(ConstantCatch)>1 & length(ConstantF)>1){stop("Cannot have both constant F and constant catch- choose one and set the other to NA!")}

  if(InitDepl > 1){
    InitDepl = 1
  }
  nyrs <- nyears + 1
  AgePart <- AgeMat + 1 # Age at first parturition = age at maturity +1 (~gestation period)

  Neq <- Ninit <- vector(length=(nages+1))
  N <- C <-  matrix(0,nrow=nyrs,ncol=(nages+1))

  Tot1P <- rep(0,length=nyrs)
  Nrep <- rep(0,length=nyrs) # number of reproductive individuals

  f0 <- (1-S1plus)/(S0*(S1plus)^(AgeMat-1))  # analytical soln for f0
  fmax <- getFecmax2(lambdaMax = lambdaMax,S1plus = S1plus,S0 = S0,AgeMat = AgeMat)

  # Equilibrium conditions (need outside of if() statement to get R0)
  Neq[1] = 1        # Age 0
  Neq[2] = S0       # Age 1
  for (a in 3:nages) {
    Neq[a] = Neq[a-1]*S1plus} #Age 2+
  Neq[nages+1]= (S0*S1plus^(nages-1))/(1-S1plus) # plus group
  R0 <- K1plus/sum(Neq[2:(nages+1)])    # numerical soln


  # Initial conditions, equilibrium
  if(InitDepl == 1){ # pop starts at equilibrium
    N[1,] <- Neq * R0
    Tot1P[1] <- K1plus
    Nrep[1] <- sum(N[1,AgePart:(nages+1)])

    # Initial conditions, non-equilibrium
  } else{             # pop starts at InitDepl*K
    NPROut <- NPR(S0 = S0,S1plus = S1plus,nages = nages, AgeMat = AgeMat, f=0)
    (unpr <- NPROut$npr)
    P1r <- NPROut$P1r # 1+ nums per recruit

    E = getF(f.start = 0.5,
             S0.w = S0,
             S1plus.w = S1plus,
             nages.w = nages,
             K1plus.w = K1plus,
             AgeMat.w = AgeMat,
             InitDepl.w = InitDepl,
             z.w = z,
             lambdaMax.w = lambdaMax,
             unpr.w = unpr)

    Ninit[1] <- 1 # N_exploited; Age 0
    Ninit[2] <- S0 # Age 1
    for(a in 3:nages){
      Ninit[a] <- S0*(S1plus*(1-E))^(a-2)
    }
    Ninit[nages+1] <- (S0*(S1plus*(1-E))^(nages-1))/(1-(S1plus*(1-E)))

    #-----
    Fec0 <- 1.0/unpr
    A <- (fmax - Fec0) / Fec0

    RF <- getRF(E = E,S0 = S0,S1plus = S1plus,nages = nages,K1plus = K1plus,AgeMat = AgeMat,z = z,A = A,unpr = unpr)
    InitNumsAtAge <- Ninit*RF    # Initial nums at age
    PropsAtAge <- InitNumsAtAge/sum(InitNumsAtAge) # Proportions at age

    n0 <- PropsAtAge[1]/sum(PropsAtAge[-1]) # Proportion of the pop that is age 0
    N0.fished <- InitDepl*K1plus*n0 # Number of age 0 individuals @ the start

    N[1,] <- N0.fished*Ninit

    Tot1P[1] <- sum(N[1,2:(nages+1)])
    Nrep[1] <- sum(N[1,AgePart:(nages+1)])
  } #end initial conditions

  if(length(ConstantCatch)>1){
    for (Yr in 1:nyears){
      MortE <- min(ConstantCatch[Yr]/Tot1P[Yr],0.99) #bycatch mortality rate

      N[Yr+1,2] <- N[Yr,1]*S0
      N[Yr+1,3:(nages+1)] <- N[Yr,2:nages]*(1-MortE)*S1plus
      N[Yr+1,(nages+1)] <- (N[Yr,nages]+N[Yr,nages+1])*(1-MortE)*S1plus
      Tot1P[Yr+1] <- sum(N[Yr+1,2:(nages+1)])
      Nrep[Yr+1] <- sum(N[Yr+1,AgePart:(nages+1)])
      N[Yr+1,1] <- Nrep[Yr+1]*(f0+(fmax-f0)*(1-(Tot1P[Yr+1]/K1plus)^z))
    }
  }else{
    sel=1
    for (Yr in 1:nyears){
      MortE <- ConstantF[Yr]*sel

      N[Yr+1,2] <- N[Yr,1]*S0
      N[Yr+1,3:(nages+1)] <- N[Yr,2:nages]*(1-MortE)*S1plus
      N[Yr+1,(nages+1)] <- (N[Yr,nages]+N[Yr,nages+1])*(1-MortE)*S1plus
      Tot1P[Yr+1] <- sum(N[Yr+1,2:(nages+1)])
      Nrep[Yr+1] <- sum(N[Yr+1,AgePart:(nages+1)])
      N[Yr+1,1] <- Nrep[Yr+1]*(f0+(fmax-f0)*(1-(Tot1P[Yr+1]/K1plus)^z)) # rec
    }
  }
  N <- N[-nyrs,]
  Tot1P <- Tot1P[-nyrs]
  return(list(TotalPop=Tot1P,N=N))
}
