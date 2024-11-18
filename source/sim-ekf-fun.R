
nl.ekf.try <- function(psi,
                       gamma,
                       times,
                       V,
                       E,
                       y,
                       hmax){
  nl.res <- try(quiet(get.nl.ekf(psi = psi,
                                 gamma = gamma,
                                 times = times,
                                 V = V,
                                 E = E,
                                 y = y,
                                 hmax = hmax)), silent = TRUE)
  if(inherits(nl.res,'try-error') |
     is.nan(nl.res) |
     is.infinite(nl.res)){
    return(1e8)
  }
  return(nl.res)
}

gnl.ekf.try <- function(psi,
                        gamma,
                        times,
                        V,
                        E,
                        y,
                        hmax){
  gnl.res <- try(quiet(get.gnl.ekf(psi = psi,
                                   gamma = gamma,
                                   times = times,
                                   V = V,
                                   E = E,
                                   y = y,
                                   hmax = hmax)), silent = TRUE)
  if(inherits(gnl.res,'try-error') |
     sum(is.nan(gnl.res)) > 0 |
     sum(is.infinite(gnl.res)) > 0){
    return(rep(1e8, length(psi)))
  }
  return(gnl.res)
}





get.nl.ekf <- function(psi,
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
  
  nTP <- length(times)
  
  xP <- matrix(data = NA,
               nrow = nTP, # nStates*(nStates - 1)/2 + 2*nStates = 14
               ncol = 14)
  
  xP[1,] <- get.xP.steady(t = max(times),
                          parms = c(theta, gamma))
  
  mu <- rep(NA, nTP)
  s <- rep(NA, nTP)
  mu[1] <- xP[1,2]*gs*nu*(V[1] - E)
  s[1] <- s2
  
  # Dt <- diff(times)[1]
  Dt <- .1
  
  for(k in 2:nTP){
    # cat(paste0("t = ", times[k], "ms\n"))
    xP[k,] <- try(vode(y = xP[k-1,], # nStates*(nStates - 1)/2 + nStates = 10
                       times = times[c(k-1,k)],
                       func = get.dxdP,
                       jacfunc = get.Jf.xP,
                       mf = -21,
                       hmax = diff(times[c(k-1,k)]), # max(tps), diff(times[c(k-1,k)])
                       parms = c(theta, gamma),
                       rtol = 1e-6),
                  silent = TRUE)[2,-1]
    
    mu[k] <- xP[k,2]*gs*nu*(V[k] - E) # "O" = 2
    s[k] <- gs^2*nu*(V[k] - E)^2*xP[k,9]*Dt + s2
    
    xP[k, 1:4] <- xP[k, 1:4] + (gs*(V[k] - E)*xP[k, c(6, 9, 10, 11)]*Dt/s[k]) * (y[k] - mu[k])
    xP[k, 5:14] <- xP[k, 5:14] - gs^2*(V[k] - E)^2*Dt^2/s[k]*xP[k, c(6,6,6,6, 
                                                                     9,9,9, 
                                                                     10,10, 
                                                                     11)] * xP[k, c(6, 9, 10, 11,
                                                                                    9, 10, 11,
                                                                                    10, 11,
                                                                                    11)]
  }
  
  nl <- as.numeric(sum(log(s)) +  t(y - mu) %*% ((y - mu)/s))
  return(nl)
}


get.gnl.ekf <- function(psi,
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
  
  nTP <- length(times)
  
  xP.wg <- matrix(data = NA,
                  nrow = nTP,
                  ncol = 126)
  
  mu <- rep(NA, nTP)
  s <- rep(NA, nTP)
  
  dx_th.idx <- tail(which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) %in% 1:4), 32)
  dP_th.idx <- tail(which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) %in% c(0,5:14)), 80)
  # xP.wg[1,1:4] <- get.x.steady(t = max(times),
  #                              parms = c(theta, gamma))
  # xP.wg[1,5:14] <- rep(0, 10)
  # xP.wg[1,dx_th.idx] <- get.x.steady_djs(t = max(times),
  #                                        parms = c(theta, gamma))
  # xP.wg[1,dP_th.idx] <- rep(0, 80)
  
  xP.wg[1,] <- c(get.xP.steady(t = max(times),
                               parms = c(theta, gamma)),
                 get.xP.steady_djs(t = max(times),
                                   parms = c(theta, gamma)))
  
  
  mu[1] <- xP.wg[1,2]*gs*nu*(V[1] - E)
  s[1] <- s2
  
  dmu_theta <- ds_theta <- matrix(data = NA,
                                  nrow = nTP,
                                  ncol = length(theta))
  dmu_gs <- dmu_nu <- ds_gs <- ds_s2 <- ds_nu <- rep(NA, nTP)
  dmu_s2 <- rep(0, nTP)
  
  dmu_theta[1, ] <- xP.wg[1, which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) == 2)[-1]]*gs*nu*(V[1] - E)
  dmu_gs[1] <- xP.wg[1,2]*nu*(V[1] - E)
  dmu_nu[1] <- xP.wg[1,2]*gs*(V[1] - E)
  
  ds_theta[1,] <- xP.wg[1, which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) == 9)[-1]]*gs^2*nu*(V[1] - E)^2*diff(times)[1]
  ds_gs[1] <- 2*gs*nu*(V[1] - E)^2*xP.wg[1,9]*diff(times)[1]
  ds_s2[1] <- 1
  ds_nu[1] <- gs^2*(V[1] - E)^2*xP.wg[1,9]*diff(times)[1]
  
  # Dt <- diff(times)[1]
  Dt <- .1
  
  for(k in 2:nTP){
    
    cat(paste0("t = ", times[k], "ms\n"))
    
    xP.wg[k,] <- try(vode(y = xP.wg[k-1,], # nStates*(nStates - 1)/2 + nStates = 10
                          times = times[c(k-1,k)],
                          func = get.dxdP.wg,
                          jacfunc = get.Jf.xP.wg,
                          mf = -21,
                          hmax = diff(times[c(k-1,k)]), # max(tps),
                          parms = c(theta, gamma),
                          rtol = c(rep(1e-6, 4),
                                   rep(1e-6, 10),
                                   rep(1e-6, 112))),
                     silent = TRUE)[2,-1]
    
    mu[k] <- xP.wg[k,2]*gs*nu*(V[k] - E)
    s[k] <- gs^2*nu*(V[k] - E)^2*xP.wg[k,9]*Dt + s2
    
    dmu_theta[k, ] <- xP.wg[k, which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) == 2)[-1]]*gs*nu*(V[k] - E)
    dmu_gs[k] <- xP.wg[k,2]*nu*(V[k] - E)
    dmu_nu[k] <- xP.wg[k,2]*gs*(V[k] - E)
    
    ds_theta[k,] <- xP.wg[k, which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) == 9)[-1]]*gs^2*nu*(V[k] - E)^2*Dt
    ds_gs[k] <- 2*gs*nu*(V[k] - E)^2*xP.wg[k,9]*Dt
    ds_s2[k] <- 1
    ds_nu[k] <- gs^2*(V[k] - E)^2*xP.wg[k,9]*Dt
    
    
    K <- gs*(V[k] - E)*xP.wg[k, c(6, 9, 10, 11)]*Dt/s[k]
    
    
    # dK_theta <- gs*(V[k] - E)*Dt*(sapply(c(6, 9, 10, 11), 
    #                                                  FUN = function(i){
    #                                                    xP.wg[k, which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) == i)[-1]]
    #                                                  })/s[k]
    #                                           - xP.wg[k, c(6, 9, 10, 11)]/(s[k]^2)*ds_theta[k,])
    dK_theta <- gs*(V[k] - E)*Dt*(matrix(data = c(xP.wg[k, c(20,  34,  48,  62,  76,  90, 104, 118,
                                                             23,  37,  51,  65,  79,  93, 107, 121,
                                                             24,  38,  52,  66,  80,  94, 108, 122,
                                                             25,  39,  53,  67,  81,  95, 109, 123)]), # which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) == 11)[-1]
                                         ncol = 4, byrow = FALSE)/s[k]
                                  - xP.wg[k, c(6, 9, 10, 11)]/(s[k]^2)*ds_theta[k,])
    
    
    dK_gs <- (V[k] - E)*xP.wg[k, c(6, 9, 10, 11)]*Dt*(1/s[k]  -gs/(s[k]^2)*ds_gs[k])
    dK_nu <- -gs*(V[k] - E)*xP.wg[k, c(6, 9, 10, 11)]*Dt/(s[k]^2)*ds_nu[k]
    dK_s2 <- -gs*(V[k] - E)*xP.wg[k, c(6, 9, 10, 11)]*Dt/(s[k]^2)*ds_s2[k]
    
    
    xP.wg[k, dx_th.idx] <- xP.wg[k, dx_th.idx] + c(t(dK_theta) * (y[k] - mu[k])) - c(K %*% t(dmu_theta[k,])) # matrix(K, nrow = 1) %*% matrix(y[k] - dmu_theta[k,], ncol = 1)
    
    
    # xP.wg[k, dP_th.idx] <- (xP.wg[k, dP_th.idx] 
    #                         - gs^2*(V[k] - E)^2*Dt^2*c(
    #                           (
    #                             t(sapply(c(6, 9, 10, 11), 
    #                                      FUN = function(i){
    #                                        xP.wg[k, which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) == i)[-1]]
    #                                      })[,c(1,1,1,1,
    #                                            2,2,2,
    #                                            3,3,
    #                                            4)]) * xP.wg[k, c(6, 9, 10, 11,
    #                                                              9, 10, 11,
    #                                                              10, 11,
    #                                                              11)] +
    #                               t(sapply(c(6, 9, 10, 11), 
    #                                        FUN = function(i){
    #                                          xP.wg[k, which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) == i)[-1]]
    #                                        })[,c(1,2,3,4,
    #                                              2,3,4,
    #                                              3,4,
    #                                              4)]) * xP.wg[k, c(6,6,6,6, 
    #                                                                9,9,9, 
    #                                                                10,10, 
    #                                                                11)])/s[k] - 
    #                             xP.wg[k, c(6,6,6,6, 
    #                                        9,9,9, 
    #                                        10,10, 
    #                                        11)] * xP.wg[k, c(6, 9, 10, 11,
    #                                                          9, 10, 11,
    #                                                          10, 11,
    #                                                          11)] %*% t(ds_theta[k,]) /(s[k]^2)) )
    xP.wg[k, dP_th.idx] <- (xP.wg[k, dP_th.idx] 
                            - gs^2*(V[k] - E)^2*Dt^2*c(
                              (
                                matrix(data = c(xP.wg[k, c(20,  34,  48,  62,  76,  90, 104, 118,
                                                           23,  37,  51,  65,  79,  93, 107, 121,
                                                           24,  38,  52,  66,  80,  94, 108, 122,
                                                           25,  39,  53,  67,  81,  95, 109, 123)]), # which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) == 11)[-1]
                                       ncol = 8, byrow = TRUE)[c(1,1,1,1,
                                                                 2,2,2,
                                                                 3,3,
                                                                 4),] * xP.wg[k, c(6, 9, 10, 11,
                                                                                   9, 10, 11,
                                                                                   10, 11,
                                                                                   11)] +
                                  matrix(data = c(xP.wg[k, c(20,  34,  48,  62,  76,  90, 104, 118,
                                                             23,  37,  51,  65,  79,  93, 107, 121,
                                                             24,  38,  52,  66,  80,  94, 108, 122,
                                                             25,  39,  53,  67,  81,  95, 109, 123)]), # which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) == 11)[-1]
                                         ncol = 8, byrow = TRUE)[c(1,2,3,4,
                                                                   2,3,4,
                                                                   3,4,
                                                                   4),] * xP.wg[k, c(6,6,6,6, 
                                                                                     9,9,9, 
                                                                                     10,10, 
                                                                                     11)])/s[k] - 
                                xP.wg[k, c(6,6,6,6, 
                                           9,9,9, 
                                           10,10, 
                                           11)] * xP.wg[k, c(6, 9, 10, 11,
                                                             9, 10, 11,
                                                             10, 11,
                                                             11)] %*% t(ds_theta[k,]) /(s[k]^2)) )
    
    
    xP.wg[k, 1:4] <- xP.wg[k, 1:4] + K * (y[k] - mu[k])
    xP.wg[k, 5:14] <- xP.wg[k, 5:14] - gs^2*(V[k] - E)^2*Dt^2/s[k]*xP.wg[k, c(6,6,6,6, 
                                                                              9,9,9, 
                                                                              10,10, 
                                                                              11)] * xP.wg[k, c(6, 9, 10, 11,
                                                                                                9, 10, 11,
                                                                                                10, 11,
                                                                                                11)]
  }
  
  gnl_theta <- as.numeric(colSums(ds_theta/s)
                          -as.numeric(2*t(dmu_theta) %*% ((y - mu)/s))
                          -as.numeric(t(y - mu) %*% Diagonal(length(s), 1/s^2) %*% (ds_theta*(y - mu))))*theta
  
  gnl_gs <- as.numeric(sum(ds_gs/s)
                       -as.numeric(2*t(dmu_gs) %*% ((y - mu)/s))
                       -as.numeric(t(y - mu) %*% Diagonal(length(s), 1/s^2) %*% (ds_gs*(y - mu))))*gs
  
  gnl_s2 <- as.numeric(sum(1/s) - t(y - mu) %*% ((y - mu)/s^2))*s2
  gnl_nu <- as.numeric(sum(ds_nu/s)
                       -as.numeric(2*t(dmu_nu) %*% ((y - mu)/s))
                       -as.numeric(t(y - mu) %*% Diagonal(length(s), 1/s^2) %*% (ds_nu*(y - mu))))*exp(tail(psi,1))
  
  return(c(gnl_theta,
           gnl_gs,
           gnl_s2,
           gnl_nu))
  
}


get.pred.ekf <- function (psi, gamma, times, V, E, y, hmax) {
  phi <- exp(psi)
  theta <- head(phi, -3)
  gs <- phi[length(theta) + 1]
  s2 <- phi[length(theta) + 2]
  nu <- 1 + phi[length(theta) + 3]
  nTP <- length(times)
  xP <- matrix(data = NA, nrow = nTP, ncol = 14)
  xP[1, ] <- c(get.x.steady(t = max(times), parms = c(theta, 
                                                      gamma)), rep(0, 10))
  mu <- rep(NA, nTP)
  s <- rep(NA, nTP)
  mu[1] <- xP[1, 2] * gs * nu * (V[1] - E)
  s[1] <- s2
  Dt <- diff(times)[1]
  for (k in 2:nTP) {
    xP[k, ] <- try(vode(y = xP[k - 1, ], 
                        times = times[c(k - 1, k)], 
                        func = get.dxdP, 
                        jacfunc = get.Jf.xP, 
                        mf = -21, 
                        hmax = diff(times[c(k - 1, k)]), parms = c(theta, 
                                                                   gamma), rtol = 1e-06), silent = TRUE)[2, -1]
    mu[k] <- xP[k, 2] * gs * nu * (V[k] - E)
    s[k] <- gs^2 * nu * (V[k] - E)^2 * xP[k, 9] * Dt + s2
    xP[k, 1:4] <- xP[k, 1:4] + (gs * (V[k] - E) * xP[k, c(6, 
                                                          9, 10, 11)] * Dt/s[k]) * (y[k] - mu[k])
    xP[k, 5:14] <- xP[k, 5:14] - gs^2 * (V[k] - E)^2 * Dt^2/s[k] * 
      xP[k, c(6, 6, 6, 6, 9, 9, 9, 10, 10, 11)] * xP[k, 
                                                     c(6, 9, 10, 11, 9, 10, 11, 10, 11, 11)]
  }
  
  return(cbind(mu,s))
}


