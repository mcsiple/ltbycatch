#' Generate one marine mammal population trajectory
#'
#' @description Generates one trajectory for a marine mammal population, starting at InitDepl
#' @param S0 calf survival rate
#' @param S1plus survival rate for animals 1 and older
#' @param K1plus pre-exploitation population size (1+ component of pop)
#' @param AgeMat age at maturity
#' @param InitDepl starting depletion level
#' @param ConstantCatch vector(length=nyears) total bycatch each year
#' @param ConstantF vector(length=nyears) rate of bycatch each year
#' @param z compensation (2.39 for bowheads, set as default for others)
#' @param nyears number of years to simulate
#' @param nages "max" age, treated as plus group
#' @param lambdaMax max intrinsic growth rate
#' @return a list containing one matrix N and one vector TotalPop
#'  N has dimensions nyears (rows) x nages (columns) and contains the number of individuals at each age. TotalPop, a vector of length nyears, contains all the 1+ individuals in the population.

# Note, nages = PlusGroupAge, and PlusGroupAge can = AgeMat+2 without losing accuracy (per AEP 11/30/18)

Dynamics <- function(S0, S1plus, K1plus, AgeMat, InitDepl, ConstantCatch=NA, ConstantF=NA, z, nyears, nages, lambdaMax){
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
