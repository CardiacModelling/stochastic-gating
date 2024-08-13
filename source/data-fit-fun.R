reversal_potential <- function(temperature){
      # Calculates the reversal potential for Potassium ions, using Nernst's
    # equation for a given ``temperature`` in degrees celsius and the internal
    # and external [K]+ concentrations used in the experiments.

    T = 273.15 + temperature
    F = 96485
    R = 8314
    K_i = 130
    k_o = 4
    return(((R*T)/F) * log(k_o/K_i))
}

temperature <- function(cell){
  # Returns the temperature (in degrees Celsius) for the given integer index
  # ``cell``.
  temperatures = c( 21.3,    # 16713003
                    21.4,    # 16715049
                    21.8,    # 16708016
                    21.7,    # 16708060
                    21.4,    # 16713110
                    21.7,    # 16708118
                    21.2,    # 16704007
                    21.6,    # 16704047
                    21.4)   # 16707014
  
  return(temperatures[cell])
}

##############################################
##### FUNCTIONS FOR Michaelis_Menten MODEL
##############################################

mm_nl <- function(psi,
                  lS,
                  K0 # mmol/L [K]o EXTRACELLULAR POTASSIUM CONCENTRATION
){
  
  lS_hat <- log(psi[1]) - log(1 + psi[2]/K0)
  
  return(as.numeric(-2*t(lS)%*%lS_hat + t(lS_hat)%*%lS_hat))
}


mm_gnl <- function(psi,
                   lS,
                   K0 # mmol/L [K]o EXTRACELLULAR POTASSIUM CONCENTRATION
){
  
  lS_hat <- log(psi[1]) - log(1 + psi[2]/K0)
  
  glS_hat_1 <- rep(1/psi[1], length(K0))
  glS_hat_2 <- -1/(K0*(1 + psi[2]/K0))
  glS_hat <- cbind(glS_hat_1, glS_hat_2)
  
  
  mmgnl <- as.numeric(-2*t(lS)%*%glS_hat)+ as.numeric(2*t(glS_hat)%*%lS_hat)
  
  return(mmgnl)
}


mm_fit <- function(psi,
                   K0 # mmol/L [K]o EXTRACELLULAR POTASSIUM CONCENTRATION
){
  
  lS_hat <-  log(psi[1]) - log(1 + psi[2]/K0)
  
  return(lS_hat)
}


mm_s2 <- function(psi,
                  lS,
                  K0 # mmol/L [K]o EXTRACELLULAR POTASSIUM CONCENTRATION
){
  
  lS_hat <- log(psi[1]) - log(1 + psi[2]/K0)
  
  s2 <- t(lS - lS_hat) %*% (lS - lS_hat)/(length(K0) - 1)
  
  return(as.numeric(s2))
}

