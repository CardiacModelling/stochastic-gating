get.nl <- function(psi,
                   gamma,
                   times,
                   V,
                   E,
                   y,
                   hmax){

  phi <- exp(psi)

  theta <- head(phi, -3)
  gs <- phi[length(theta) + 1]
  s2 <- phi[length(theta) + 2]
  nu <- 1 + phi[length(theta) + 3]
  
  Dt <- .1

  xP <- try(vode(y = get.xP.steady(t = max(times),
                                    parms = c(theta, gamma)), # nStates*(nStates - 1)/2 + nStates = 10
                 times = times,
                 func = get.dxdP,
                 jacfunc = get.Jf.xP,
                 mf = -21,
                 hmax = hmax, # max(tps),
                 parms = c(theta, gamma),
                 rtol = 1e-6),
            silent = TRUE)[,-1]

  mu <- xP[,"O"]*gs*nu*(V - E)
  s <- gs^2*nu*(V - E)^2*xP[,9]*Dt + s2

  nl <- as.numeric(sum(log(s)) +  t(y - mu) %*% ((y - mu)/s))
  return(nl)
}


get.gnl <- function(psi,
                    gamma,
                    times,
                    V,
                    E,
                    y,
                    hmax){

  phi <- exp(psi)

  theta <- head(phi, -3)
  gs <- phi[length(theta) + 1]
  s2 <- phi[length(theta) + 2]
  nu <- 1 + phi[length(theta) + 3]
  
  Dt <- .1

  xP.wg <- try(vode(y = c(get.xP.steady(t = max(times),
                                    parms = c(theta, gamma)),
                          get.xP.steady_djs(t = max(times),
                                           parms = c(theta, gamma))),
                    times = times,
                    func = get.dxdP.wg,
                    jacfunc = get.Jf.xP.wg,
                    mf = -21,
                    hmax = hmax, # max(tps),
                    parms = c(theta, gamma),
                    rtol = c(rep(1e-6, 4),
                             rep(1e-6, 10),
                             rep(1e-6, 112))),
               silent = TRUE)[,-1]

  mu <- xP.wg[,"O"]*gs*nu*(V - E)
  s <- gs^2*nu*(V - E)^2*xP.wg[,9]*Dt + s2

  dmu_theta <- xP.wg[, which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) == 2)[-1]]*gs*nu*(V - E)
  dmu_gs <- xP.wg[,"O"]*nu*(V - E)
  dmu_nu <- xP.wg[,"O"]*gs*(V - E)

  ds_theta <- xP.wg[, which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) == 9)[-1]]*gs^2*nu*(V - E)^2*Dt
  ds_gs <- 2*gs*nu*(V - E)^2*xP.wg[,9]*Dt
  ds_s2 <- 1
  ds_dnu <- gs^2*(V - E)^2*xP.wg[,9]*Dt

  gnl_theta <- as.numeric(colSums(ds_theta/s)
                          -as.numeric(2*t(dmu_theta) %*% ((y - mu)/s))
                          -as.numeric(t(y - mu) %*% Diagonal(length(s), 1/s^2) %*% (ds_theta*(y - mu))))*theta

  gnl_gs <- as.numeric(sum(ds_gs/s)
                       -as.numeric(2*t(dmu_gs) %*% ((y - mu)/s))
                       -as.numeric(t(y - mu) %*% Diagonal(length(s), 1/s^2) %*% (ds_gs*(y - mu))))*gs

  gnl_s2 <- as.numeric(sum(1/s) - t(y - mu) %*% ((y - mu)/s^2))*s2
  gnl_nu <- as.numeric(sum(ds_dnu/s)
                       -as.numeric(2*t(dmu_nu) %*% ((y - mu)/s))
                       -as.numeric(t(y - mu) %*% Diagonal(length(s), 1/s^2) %*% (ds_dnu*(y - mu))))*exp(tail(psi,1))

  return(c(gnl_theta,
           gnl_gs,
           gnl_s2,
           gnl_nu))

}

# get.nl <- function(psi,
#                    gamma,
#                    times,
#                    V,
#                    E,
#                    y,
#                    hmax){
# 
#   phi <- exp(psi)
# 
#   theta <- head(phi, -3)
#   gs <- phi[length(theta) + 1]
#   s2 <- phi[length(theta) + 2]
#   nu <- 1 + phi[length(theta) + 3]
# 
#   xP <- try(vode(y = c(get.x.steady(t = max(times),
#                                     parms = c(theta, gamma)),
#                        rep(0, 10)), # nStates*(nStates - 1)/2 + nStates = 10
#                  times = times,
#                  func = get.dxdP,
#                  jacfunc = get.Jf.xP,
#                  mf = -21,
#                  hmax = hmax, # max(tps),
#                  parms = c(theta, gamma),
#                  rtol = 1e-6),
#             silent = TRUE)[,-1]
# 
#   mu <- xP[,"O"]*gs*nu*(V - E)
#   # s <- gs^2*nu*(V - E)^2*xP[,9]*diff(times)[1] + s2
#   Dt <- c(0, diff(times))
#   # Dt <- .1
#   s <- gs^2*nu*(V - E)^2*xP[,9]*Dt + s2
# 
#   nl <- as.numeric(sum(log(s)) +  t(y - mu) %*% ((y - mu)/s))
#   return(nl)
# }
# 
# get.gnl <- function(psi,
#                     gamma,
#                     times,
#                     V,
#                     E,
#                     y,
#                     hmax){
# 
#   phi <- exp(psi)
# 
#   theta <- head(phi, -3)
#   gs <- phi[length(theta) + 1]
#   s2 <- phi[length(theta) + 2]
#   nu <- 1 + phi[length(theta) + 3]
# 
#   xP.wg <- try(vode(y = c(get.x.steady(t = max(times),
#                                        parms = c(theta, gamma)),
#                           rep(0, 10), # nStates*(nStates - 1)/2 + nStates
#                           get.x.steady_djs(t = max(times),
#                                            parms = c(theta, gamma)),
#                           rep(0, 80)),
#                     times = times,
#                     func = get.dxdP.wg,
#                     jacfunc = get.Jf.xP.wg,
#                     mf = -21,
#                     hmax = hmax, # max(tps),
#                     parms = c(theta, gamma),
#                     rtol = c(rep(1e-6, 4),
#                              rep(1e-6, 10),
#                              rep(1e-6, 112))),
#                silent = TRUE)[,-1]
# 
#   mu <- xP.wg[,"O"]*gs*nu*(V - E)
#   # s <- gs^2*nu*(V - E)^2*xP.wg[,9]*diff(times)[1] + s2
#   Dt <- c(0, diff(times))
#   # Dt <- .1
#   s <- gs^2*nu*(V - E)^2*xP.wg[,9]*Dt + s2
# 
#   dmu_theta <- xP.wg[, which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) == 2)[-1]]*gs*nu*(V - E)
#   dmu_gs <- xP.wg[,"O"]*nu*(V - E)
#   dmu_nu <- xP.wg[,"O"]*gs*(V - E)
# 
#   # ds_theta <- xP.wg[, which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) == 9)[-1]]*gs^2*nu*(V - E)^2*diff(times)[1]
#   # ds_gs <- 2*gs*nu*(V - E)^2*xP.wg[,9]*diff(times)[1]
#   # ds_s2 <- 1
#   # ds_dnu <- gs^2*(V - E)^2*xP.wg[,9]*diff(times)[1]
#   ds_theta <- xP.wg[, which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) == 9)[-1]]*gs^2*nu*(V - E)^2*Dt
#   ds_gs <- 2*gs*nu*(V - E)^2*xP.wg[,9]*Dt
#   ds_s2 <- 1
#   ds_dnu <- gs^2*(V - E)^2*xP.wg[,9]*Dt
# 
#   gnl_theta <- as.numeric(colSums(ds_theta/s)
#                           -as.numeric(2*t(dmu_theta) %*% ((y - mu)/s))
#                           -as.numeric(t(y - mu) %*% Diagonal(length(s), 1/s^2) %*% (ds_theta*(y - mu))))*theta
# 
#   gnl_gs <- as.numeric(sum(ds_gs/s)
#                        -as.numeric(2*t(dmu_gs) %*% ((y - mu)/s))
#                        -as.numeric(t(y - mu) %*% Diagonal(length(s), 1/s^2) %*% (ds_gs*(y - mu))))*gs
# 
#   gnl_s2 <- as.numeric(sum(1/s) - t(y - mu) %*% ((y - mu)/s^2))*s2
#   gnl_nu <- as.numeric(sum(ds_dnu/s)
#                        -as.numeric(2*t(dmu_nu) %*% ((y - mu)/s))
#                        -as.numeric(t(y - mu) %*% Diagonal(length(s), 1/s^2) %*% (ds_dnu*(y - mu))))*exp(tail(psi,1))
# 
#   return(c(gnl_theta,
#            gnl_gs,
#            gnl_s2,
#            gnl_nu))
# 
# }



get.nl.nus2 <- function(psi,
                        g,
                        gamma,
                        times,
                        V,
                        E,
                        y,
                        xP.wg){
  
  phi <- exp(psi)
  s2 <- phi[1]
  nu <- 1 + phi[2]
  
  mu <- xP.wg[,"O"]*g*(V - E)
  # s <- g^2*(V - E)^2*xP.wg[,9]/nu*diff(times)[1] + s2
  # Dt <- c(0, diff(times))
  Dt <- .1
  s <- g^2*(V - E)^2*xP.wg[,9]/nu*Dt + s2
  
  nl <- as.numeric(sum(log(s)) +  t(y - mu) %*% ((y - mu)/s))
  return(nl)
}

get.gnl.nus2 <- function(psi,
                         g,
                         gamma,
                         times,
                         V,
                         E,
                         y,
                         xP.wg){
  
  phi <- exp(psi)
  s2 <- phi[1]
  nu <- 1 + phi[2]
  
  mu <- xP.wg[,"O"]*g*(V - E)
  # s <- g^2*(V - E)^2*xP.wg[,9]/nu*diff(times)[1] + s2
  # Dt <- c(0, diff(times))
  Dt <- .1
  s <- g^2*(V - E)^2*xP.wg[,9]/nu*Dt + s2
  
  dmu_s2 <- 0
  # ds_dnu <- -(g^2*(V - E)^2*xP.wg[,9]/nu^2)*diff(times)[1]
  ds_dnu <- -(g^2*(V - E)^2*xP.wg[,9]/nu^2)*Dt
  
  gnl_s2 <- as.numeric(sum(1/s) - t(y - mu) %*% ((y - mu)/s^2))*s2
  gnl_nu <- as.numeric(sum(ds_dnu/s)
                       -as.numeric(t(y - mu) %*% Diagonal(length(s), 1/s^2) %*% (ds_dnu*(y - mu))))*exp(tail(psi,1))
  
  return(c(gnl_s2,
           gnl_nu))
}


# get.pred <- function(psi,
#                      gamma,
#                      times,
#                      V,
#                      E,
#                      hmax){
# 
#   phi <- exp(psi)
# 
#   theta <- head(phi, -3)
#   gs <- phi[length(theta) + 1]
#   s2 <- phi[length(theta) + 2]
#   nu <- 1 + phi[length(theta) + 3]
# 
#   xP <- try(vode(y = c(get.x.steady(t = max(times),
#                                     parms = c(theta, gamma)),
#                        rep(0, 10)), # nStates*(nStates - 1)/2 + nStates = 10
#                  times = times,
#                  func = get.dxdP,
#                  jacfunc = get.Jf.xP,
#                  mf = -21,
#                  hmax = hmax, # max(tps),
#                  parms = c(theta, gamma),
#                  rtol = 1e-6),
#             silent = TRUE)[,-1]
#   
#   # Dt <- c(0, diff(times))
#   Dt <- .1
# 
#   mu <- xP[,"O"]*gs*nu*(V - E)
#   s <- gs^2*nu*(V - E)^2*xP[,9]*Dt + s2
# 
#   return(cbind(mu, s))
# }


