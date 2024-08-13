get.pred.X <- function (psi, gamma, times, V, E, hmax){
  phi <- exp(psi)
  theta <- head(phi, -3)
  gs <- phi[length(theta) + 1]
  s2 <- phi[length(theta) + 2]
  nu <- 1 + phi[length(theta) + 3]
  xP <- try(vode(y = c(get.x.steady(t = max(times), parms = c(theta, 
                                                              gamma)), rep(0, 10)), 
                 times = times, func = get.dxdP, 
                 jacfunc = get.Jf.xP, 
                 mf = -21, 
                 hmax = hmax, 
                 parms = c(theta, gamma),
                 rtol = 1e-06), silent = TRUE)[, -1]
  return(xP)
}

