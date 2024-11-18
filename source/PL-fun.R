
FIM_parabola <- function(x, 
                         a, 
                         H,
                         gamma,
                         times,
                         V,
                         E,
                         y,
                         hmax){
  
  return(as.numeric(nl.try(psi = a,
                           gamma = gamma,
                           times = times,
                           V = V,
                           E = E,
                           y = y,
                           hmax = hmax) + .5*t(x - a)%*% H %*% (x - a)))
  
}

get.PL <- function(phi,
                   i,
                   IC_i,
                   gamma,
                   times,
                   V,
                   E,
                   y,
                   hmax){

  phi.c <- phi

  PL.values <- NULL
  PL.allres <- list()
  for (nit in 1:length(IC_i)) {

    phi_i.c <- IC_i[nit]

    cat(paste0("\tStep ", nit, "\t PL for variable ", i, "\ttheta[", i, "] = ", phi_i.c, " \n"))
    cat(paste0("\tStep ", nit, "\t PL for variable ", i, ";\tComputing Profile Likelihood...\n"))

    res.fit.PL <- try(optim(par = phi.c[-i], # psi0
                            fn = nl.PL.try,
                            gr = gnl.PL.try,
                            gamma = gamma,
                            times = times,
                            V = V,
                            E = E,
                            y = y,
                            hmax = hmax,
                            psi.idx = i,
                            psi.val = phi_i.c, # max(tps),
                            method = "L-BFGS-B",
                            control = list(lmm = 2500,
                                           factr = 0,
                                           maxit = 10000,
                                           trace = 1,
                                           REPORT = 1)),
                      silent = TRUE)

    phi.c[-i] <- res.fit.PL$par
    phi.c[i] <- phi_i.c

    PL.values <- rbind(PL.values, c(phi.c, PL = res.fit.PL$value))
    PL.allres[[nit]] <- res.fit.PL
  }


  res <- list()
  res$PL.values <- PL.values
  res$PL.allres <- PL.allres

  return(res)
}


nl.PL.try <- function(omega,
                      gamma,
                      times,
                      V,
                      E,
                      y,
                      hmax,
                      psi.idx,
                      psi.val){
  nl.res <- try(get.nl.PL(omega = omega,
                          gamma = gamma,
                          times = times,
                          V = V,
                          E = E,
                          y = y,
                          hmax = hmax,
                          psi.idx = psi.idx,
                          psi.val = psi.val), silent = TRUE)
  if(inherits(nl.res,'try-error') |
     is.nan(nl.res) |
     is.infinite(nl.res)){
    return(1e8)
  }
  return(nl.res)
}


gnl.PL.try <- function(omega,
                       gamma,
                       times,
                       V,
                       E,
                       y,
                       hmax,
                       psi.idx,
                       psi.val){
  gnl.res <- try(get.gnl.PL(omega = omega,
                            gamma = gamma,
                            times = times,
                            V = V,
                            E = E,
                            y = y,
                            hmax = hmax,
                            psi.idx = psi.idx,
                            psi.val = psi.val), silent = TRUE)
  if(inherits(gnl.res,'try-error') |
     sum(is.nan(gnl.res)) > 0 |
     sum(is.infinite(gnl.res)) > 0){
    return(rep(1e8, length(omega)))
  }
  return(gnl.res)
}

get.nl.PL <- function(omega,
                      gamma,
                      times,
                      V,
                      E,
                      y,
                      hmax,
                      psi.idx,
                      psi.val){
  
  
  psi <- rep(NA, length(omega) + length(psi.idx))
  psi[-psi.idx] <- omega
  psi[psi.idx] <- psi.val
  phi <- exp(psi)
  
  theta <- head(phi, -3)
  gs <- phi[length(theta) + 1]
  s2 <- phi[length(theta) + 2]
  nu <- 1 + phi[length(theta) + 3]
  
  xP <- try(vode(y = c(get.x.steady(t = max(times),
                                    parms = c(theta, gamma)),
                       rep(0, 10)), # nStates*(nStates - 1)/2 + nStates = 10
                 times = times,
                 func = get.dxdP,
                 jacfunc = get.Jf.xP,
                 mf = -21,
                 hmax = hmax, # max(tps),
                 parms = c(theta, gamma),
                 rtol = 1e-6),
            silent = TRUE)[,-1]
  
  mu <- xP[,"O"]*gs*nu*(V - E)
  s <- gs^2*nu*(V - E)^2*xP[,9]*diff(times)[1] + s2
  
  nl <- as.numeric(sum(log(s)) +  t(y - mu) %*% ((y - mu)/s))
  return(nl)
}


get.gnl.PL <- function(omega,
                       gamma,
                       times,
                       V,
                       E,
                       y,
                       hmax,
                       psi.idx,
                       psi.val){
  
  psi <- rep(NA, length(omega) + length(psi.idx))
  psi[-psi.idx] <- omega
  psi[psi.idx] <- psi.val
  phi <- exp(psi)
  
  theta <- head(phi, -3)
  gs <- phi[length(theta) + 1]
  s2 <- phi[length(theta) + 2]
  nu <- 1 + phi[length(theta) + 3]
  
  xP.wg <- try(vode(y = c(get.x.steady(t = max(times),
                                       parms = c(theta, gamma)),
                          rep(0, 10), # nStates*(nStates - 1)/2 + nStates
                          get.x.steady_djs(t = max(times),
                                           parms = c(theta, gamma)),
                          rep(0, 80)),
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
  s <- gs^2*nu*(V - E)^2*xP.wg[,9]*diff(times)[1] + s2
  
  dmu_theta <- xP.wg[, which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) == 2)[-1]]*gs*nu*(V - E)
  dmu_gs <- xP.wg[,"O"]*nu*(V - E)
  dmu_nu <- xP.wg[,"O"]*gs*(V - E)
  
  ds_theta <- xP.wg[, which(1:ncol(xP.wg) %% (ncol(xP.wg)/9) == 9)[-1]]*gs^2*nu*(V - E)^2*diff(times)[1]
  ds_gs <- 2*gs*nu*(V - E)^2*xP.wg[,9]*diff(times)[1]
  ds_s2 <- 1
  ds_dnu <- gs^2*(V - E)^2*xP.wg[,9]*diff(times)[1]
  
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
  
  gnl_all <- c(gnl_theta,
               gnl_gs,
               gnl_s2,
               gnl_nu)
  
  return(gnl_all[-psi.idx])
  
}

get.ellipse <- function(phi,
                        levels = c(.99, .95, .9),
                        fim,
                        N){
  
  levels <- sort(levels, decreasing = TRUE)
  
  p <- length(phi)
  cov.phi <- solve(fim)#/N
  t <- seq(0, 2*pi, length.out = 100)
  
  par(mar = c(.5,.5,.5,.5), mfrow = c(p,p))
  for (j in (p-1):1) {
    if(j >= 2){
      for (h in 1:(j - 1)) {
        plot.new()
      }
    }
    for (i in (j+1):p) {
      x0 <- phi[i]
      y0 <- phi[j]
      
      eigen_ij <- eigen(cov.phi[c(i,j), c(i,j)])
      max.v <- eigen_ij$vectors[,1]
      alpha <- atan(max.v[2]/max.v[1])
      
      par(mar = c(.5,.5,.5,.5))
      
      for (k in 1:length(levels)) {
        ab <-  2*sqrt(qchisq(p = levels[k], df = 2) * eigen_ij$values)
        a <- ab[1]
        b <- ab[2]
        xi <- a*cos(t)
        xj <- b*sin(t)
        
        if(k == 1){
          plot(xi*cos(alpha) - xj*sin(alpha) + x0,
               xi*sin(alpha) + xj*cos(alpha) + y0, 
               type = 'l', 
               lwd = 2,
               cex.axis = 1,
               cex.lab = 2,
               bty = 'n',
               yaxt = 'n',
               xaxt = ifelse(j == 1, "s", "n"),
               xlab = bquote(log~theta[.(i)]),
               ylab = bquote(log~theta[.(j)]))
          max.range.xi <- xi*cos(alpha) - xj*sin(alpha) + x0
          max.range.xj <- xi*sin(alpha) + xj*cos(alpha) + y0
          if(i == p){
            axis(side = 4)
            mtext(bquote(log~theta[.(j)]), 
                  side=4, 
                  line=3,
                  cex = 1.5)
          }
          if(j == 1){
            axis(side = 1)
            mtext(bquote(log~theta[.(i)]), 
                  side=1, 
                  line=3,
                  cex = 1.5)
          }
        }else{
          lines(xi*cos(alpha) - xj*sin(alpha) + x0,
                xi*sin(alpha) + xj*cos(alpha) + y0)
        }
      }
      points(x0, y0, pch = "*", cex = 3)
      
    }
    par(mar = c(.1,.1,.1,.1))
    plot.new()
    
  }
  
}


get.ellipse.axes.range <- function(phi,
                                   level = .9,
                                   fim,
                                   N){
  
  p <- length(phi)
  cov.phi <- solve(fim)#/N
  t <- seq(0, 2*pi, length.out = 100)
  
  ranges <- matrix(data = NA, nrow = p, ncol = 2)
  
  for (j in (p-1):1) {
    for (i in (j+1):p) {
      x0 <- phi[i]
      y0 <- phi[j]
      
      eigen_ij <- eigen(cov.phi[c(i,j), c(i,j)])
      max.v <- eigen_ij$vectors[,1]
      alpha <- atan(max.v[2]/max.v[1])
      
      ab <-  2*sqrt(qchisq(p = level, df = 2) * eigen_ij$values)
      a <- ab[1]
      b <- ab[2]
      xi <- a*cos(t)
      xj <- b*sin(t)
      
      range.xi <- range(xi*cos(alpha) - xj*sin(alpha) + x0)
      range.xj <- range(xi*sin(alpha) + xj*cos(alpha) + y0)
      
      ranges[i,] <- range.xi
      ranges[j,] <- range.xj
      
    }
  }
  
  return(ranges)
}


nl.PL_try <- function(omega,
         p9,
         y,
         times,
         V,
         E,
         hmax,
         psi.idx,
         psi.val){
  nl.res <- try(nl.PL(omega = omega,
                   p9 = p9,
                   y = y,
                   times = times,
                   V = V,
                   E = E,
                   hmax = hmax,
                   psi.idx = psi.idx,
                   psi.val = psi.val), silent = TRUE)
  if(inherits(nl.res,'try-error') |
     is.nan(nl.res) |
     is.infinite(nl.res)){
    return(1e8)
  }
  return(nl.res)
}

nl.PL <- function(omega,
         p9,
         y,
         times,
         V,
         E,
         hmax,
         psi.idx,
         psi.val){
  
  psi <- rep(NA, length(omega) + 1)
  psi[-psi.idx] <- omega
  psi[psi.idx] <- psi.val
  phi <- exp(psi)
  theta <- head(phi, -2)
  g <- tail(phi, 2)[1]
  s2 <- tail(phi, 1)
  
  x <- try(vode(y = get.x.steady(t = max(times),
                                 parms = c(theta,
                                           p9)),
                times = times,
                func = get.dx,
                jacfunc = get.Jf,
                mf = -21,
                hmax = hmax,
                parms = c(theta,
                          p9)),
           silent = TRUE)[,-1]
  mu <- x[,'O']*g*(V - E)
  n <- length(y)
  return(as.numeric(n*log(s2) +  t(y - mu) %*% ((y - mu)/s2)))
}


gnl.PL_try <- function(omega,
         p9,
         y,
         times,
         V,
         E,
         hmax,
         psi.idx,
         psi.val){
  gnl.res <- try(gnl.PL(omega = omega,
                     p9 = p9,
                     y = y,
                     times = times,
                     V = V,
                     E = E,
                     hmax = hmax,
                     psi.idx = psi.idx,
                     psi.val = psi.val), silent = TRUE)
  if(inherits(gnl.res,'try-error') |
     sum(is.nan(gnl.res)) > 0 |
     sum(is.infinite(gnl.res)) > 0){
    return(rep(1e8, length(omega)))
  }
  return(gnl.res)
}

gnl.PL <- function(omega,
         p9,
         y,
         times,
         V,
         E,
         hmax,
         psi.idx,
         psi.val){
  
  psi <- rep(NA, length(omega) + 1)
  psi[-psi.idx] <- omega
  psi[psi.idx] <- psi.val
  phi <- exp(psi)
  theta <- head(phi, -2)
  g <- tail(phi, 2)[1]
  s2 <- tail(phi, 1)
  
  x <- try(vode(y = c(get.x.steady(t = max(times),
                                   parms = c(theta,
                                             p9)),
                      get.x.steady_djs(t = max(times),
                                       parms = c(theta,
                                                 p9))),
                times = times,
                func = get.dx.wg,
                jacfunc = get.Jf.wg,
                mf = -21,
                hmax = hmax,
                parms = c(theta,
                          p9)),
           silent = TRUE)[,-1]
  
  mu <- x[,'O']*g*(V - E)
  dmu_th <- x[,which(1:ncol(x) %% (ncol(x)/(length(theta) + 1)) == which(colnames(x) == 'O'))[-1]]*g*(V - E)
  dmu_g <- x[,'O']*(V - E)
  
  n <- length(y)
  nl_th <- (-2*t(dmu_th) %*% ((y - mu)/s2))*theta
  nl_g <- (-2*t(dmu_g) %*% ((y - mu)/s2))*g
  nl_s2 <- (n - t(y - mu) %*% ((y - mu)/s2))
  
  return(c(nl_th, nl_g, nl_s2)[-psi.idx])
}

nl.nus2.PL.try <- function(omega,
                           g,
                           gamma,
                           times,
                           V,
                           E,
                           y,
                           xP.wg,
                           psi.idx,
                           psi.val){
  nl.res <- try(get.nl.nus2.PL(omega = omega,
                               g = g,
                               gamma = gamma,
                               times = times,
                               V = V,
                               E = E,
                               y = y,
                               xP.wg = xP.wg,
                               psi.idx = psi.idx,
                               psi.val = psi.val), silent = TRUE)
  if(inherits(nl.res,'try-error') |
     is.nan(nl.res) |
     is.infinite(nl.res)){
    return(1e8)
  }
  return(nl.res)
}

get.nl.nus2.PL <- function(omega,
                           g,
                           gamma,
                           times,
                           V,
                           E,
                           y,
                           xP.wg,
                           psi.idx,
                           psi.val){
  
  psi <- rep(NA, length(omega) + 1)
  psi[-psi.idx] <- omega
  psi[psi.idx] <- psi.val
  
  phi <- exp(psi)
  s2 <- phi[1]
  nu <- 1 + phi[2]
  
  mu <- xP.wg[,"O"]*g*(V - E)
  s <- g^2*(V - E)^2*xP.wg[,9]/nu*diff(times)[1] + s2
  
  nl <- as.numeric(sum(log(s)) +  t(y - mu) %*% ((y - mu)/s))
  return(nl)
}


gnl.nus2.PL.try <- function(omega,
                            g,
                            gamma,
                            times,
                            V,
                            E,
                            y,
                            xP.wg,
                            psi.idx,
                            psi.val){
  gnl.res <- try(get.gnl.nus2.PL(omega = omega,
                                 g = g,
                                 gamma = gamma,
                                 times = times,
                                 V = V,
                                 E = E,
                                 y = y,
                                 xP.wg = xP.wg,
                                 psi.idx = psi.idx,
                                 psi.val = psi.val), silent = TRUE)
  if(inherits(gnl.res,'try-error') |
     sum(is.nan(gnl.res)) > 0 |
     sum(is.infinite(gnl.res)) > 0){
    return(rep(1e8, length(psi)))
  }
  return(gnl.res)
}

get.gnl.nus2.PL <- function(omega,
                            g,
                            gamma,
                            times,
                            V,
                            E,
                            y,
                            xP.wg,
                            psi.idx,
                            psi.val){
  
  psi <- rep(NA, length(omega) + 1)
  psi[-psi.idx] <- omega
  psi[psi.idx] <- psi.val
  
  phi <- exp(psi)
  s2 <- phi[1]
  nu <- 1 + phi[2]
  
  
  mu <- xP.wg[,"O"]*g*(V - E)
  s <- g^2*(V - E)^2*xP.wg[,9]/nu*diff(times)[1] + s2
  
  dmu_s2 <- 0
  ds_dnu <- -(g^2*(V - E)^2*xP.wg[,9]/nu^2)*diff(times)[1]
  
  gnl_s2 <- as.numeric(sum(1/s) - t(y - mu) %*% ((y - mu)/s^2))*s2
  gnl_nu <- as.numeric(sum(ds_dnu/s)
                       -as.numeric(t(y - mu) %*% Diagonal(length(s), 1/s^2) %*% (ds_dnu*(y - mu))))*exp(tail(psi,1))
  
  return(c(gnl_s2,
           gnl_nu)[-psi.idx])
}


