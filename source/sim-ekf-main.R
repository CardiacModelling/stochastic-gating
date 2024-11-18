
cat("\nInstall/load packages\n")
inst.pkgs <- installed.packages()
l.pkgs <- c("stringr",
            "pracma",
            "Deriv",
            "deSolve",
            "Matrix",
            "mvtnorm",
            "tmvtnorm",
            "parallel",
            "nloptr",
            "spsUtil")
lapply(l.pkgs, function(pkg){
  if(!(pkg %in% rownames(inst.pkgs))){
    install.packages(pkg)
  }
})
lapply(l.pkgs, function(pkg){
  library(pkg,character.only=TRUE)
})

source("./source/ODE-fun.R")
source("./source/SDE-fun.R")
source("./source/sim-ekf-fun.R")
source("./source/sim-sampFreq-fun.R")

s2_nu.true_nSim_nCores <- commandArgs(trailingOnly = TRUE) # s2_nu.true_nSim_nCores <- "0.00001-1000-100-20"
s2.true <- as.numeric(unlist(strsplit(s2_nu.true_nSim_nCores, split = "-"))[1])
nu.true <- as.numeric(unlist(strsplit(s2_nu.true_nSim_nCores, split = "-"))[2])
nSim <- as.numeric(unlist(strsplit(s2_nu.true_nSim_nCores, split = "-"))[3])
nCores <- as.numeric(unlist(strsplit(s2_nu.true_nSim_nCores, split = "-"))[4])
pFreq <- as.numeric(unlist(strsplit(s2_nu.true_nSim_nCores, split = "-"))[5])

cat(paste0("\ts2 = ", s2.true, 
           "\tnu = ", nu.true, 
           "\tnSim = ", nSim, 
           "\tnCores = ", nCores, 
           "\tnpFreq = ", pFreq,
           "\n"))

##############
## currTime ##
##############

currTime <- Sys.time()
currTime <- gsub(" ", "_", currTime)
currTime <- gsub("-", "", currTime)
currTime <- gsub(":", "", currTime)
currTime <- paste(currTime, s2_nu.true_nSim_nCores, sep = "_")

#############
## folders ##
#############

currDir <- getwd()
setwd(currDir)

#############
## folders ##
#############

resPath <- paste("./results", sep = "")
ifelse(!dir.exists(file.path(resPath)), dir.create(file.path(resPath)), FALSE)
resUQ1 <- paste0(resPath, "/UQ1/")
ifelse(!dir.exists(file.path(resUQ1)), dir.create(file.path(resUQ1)), FALSE)
allFigsPath <- paste0(resPath, "allFigures/")
ifelse(!dir.exists(file.path(allFigsPath)), dir.create(file.path(allFigsPath)), FALSE)
resPath <- paste0(paste(resPath, currTime, sep = "/"), "/")
figPath <- paste0(resPath, "figures/")
printPath <- paste0(resPath, "print/")

ifelse(!dir.exists(file.path(resPath)), dir.create(file.path(resPath)), FALSE)
ifelse(!dir.exists(file.path(figPath)), dir.create(file.path(figPath)), FALSE)
ifelse(!dir.exists(file.path(printPath)), dir.create(file.path(printPath)), FALSE)

tps <- seq(from = 0,
           to = 8000,
           by = .1)
# tps <- seq(from = 0,
#            to = 800,
#            by = .1)
# tps <- seq(from = 0,
#            to = 300,
#            by = .1)

chs <- c("C", "O", "F", "I", "IC")
rcs <- c("C->O",
         "O->C",
         "C->F",
         "F->C",
         "O->I",
         "I->O",
         "F->I",
         "I->F",
         "IC->I",
         "I->IC",
         "C->IC",
         "IC->C")

cnstrs <- c("lambda\\[\\'IC->I\\'\\]=lambda\\[\\'C->O\\'\\]",
            "lambda\\[\\'I->IC\\'\\]=lambda\\[\\'O->C\\'\\]",
            "lambda\\[\\'C->IC\\'\\]=lambda\\[\\'O->I\\'\\]",
            "lambda\\[\\'IC->C\\'\\]=lambda\\[\\'I->O\\'\\]",
            "lambda\\[\\'F->C\\'\\]=lambda\\[\\'O->C\\'\\]",
            "lambda\\[\\'F->I\\'\\]=lambda\\[\\'O->I\\'\\]",
            "lambda\\[\\'C->F\\'\\]=lambda\\[\\'C->O\\'\\]",
            "lambda\\[\\'I->F\\'\\]=lambda\\[\\'I->O\\'\\]")

th.true <- c(2.232057e-04,
             7.007999e-02,
             3.405121e-05,
             5.450178e-02,
             8.708426e-02,
             8.259672e-03,
             5.397629e-03,
             3.238714e-02,
             c(1/6.7)/sum(c(1/6.7,
                            1/2.5)))

g.true <- 1.456307e-01

FLKR <- TRUE
ntheta <- 9
nStates <- length(chs) - 1
JSENS <- 1:(ntheta-1)

## SDE functions
cat(paste0("\t generating SDE functions...\n"))

generate.lambda(rct.lst = rcs,
                constr.lst = cnstrs,
                flickering = FLKR)

generate.get.Vlambda(ct.lst = chs,
                     rct.lst = rcs,
                     flickering = FLKR)

generate.get.x.steady(p = ntheta,
                      ct.lst = chs,
                      rct.lst = rcs,
                      flickering = FLKR)

generate.get.dx(ct.lst = chs,
                rct.lst = rcs,
                flickering = FLKR)

generate.get.Jf_y(ct.lst = chs,
                  rct.lst = rcs,
                  p = ntheta,
                  jsens = JSENS,
                  flickering = TRUE)

generate.bdiagmat.xp(n = nStates*(nStates - 1)/2 + 2*nStates,
                     p = length(JSENS))

generate.get.dxdP_dj(p = ntheta,
                     n = nStates*(nStates - 1)/2 + 2*nStates,
                     jsens = JSENS)

generate.get.x.steady_djs(p = ntheta,
                          jsens = JSENS,
                          ct.lst = chs,
                          rct.lst = rcs,
                          envir = .GlobalEnv,
                          flickering = FLKR)

generate.get.dxdP.wg(n = nStates*(nStates - 1)/2 + 2*nStates,
                     p = length(JSENS))

generate.get.Sigma(ct.lst = chs,
                   rct.lst = rcs)

generate.get.xP.steady(ct.lst = chs)
generate.get.xP.steady_djs(p = ntheta, 
                           jsens = JSENS)

## ODE functions
cat(paste0("\t generating ODE functions...\n"))

generate.get.dx.wg(n = nStates,
                   p = length(JSENS),
                   ct.lst = chs,
                   rct.lst = rcs,
                   envir = .GlobalEnv,
                   flickering = FLKR)

generate.get.Jf.wg(n = nStates)

generate.bdiagmat(n = nStates,
                  p = length(JSENS))

generate.get.dx_dj_di(p = ntheta,
                      jsens = JSENS,
                      n = nStates,
                      ct.lst = chs,
                      rct.lst = rcs,
                      envir = .GlobalEnv)

generate.get.dx_dj(p = ntheta,
                   jsens = JSENS,
                   n = nStates,
                   ct.lst = chs,
                   rct.lst = rcs,
                   envir = .GlobalEnv)

generate.get.dx(ct.lst = chs,
                rct.lst = rcs,
                flickering = FLKR)

generate.nl_try(p = ntheta,
                jsens = JSENS)

generate.gnl_try(p = ntheta,
                 jsens = JSENS)

generate.nl(p = ntheta,
            jsens = JSENS)

generate.gnl(p = ntheta,
             jsens = JSENS)

generate.fim(p = ntheta,
             jsens = JSENS)

cat(paste0("\t DONE.\n"))

Vsin <- sapply(tps, function(t){get.voltage(t)})
E.true <- -88.35746
gs.true <- g.true/nu.true

psi.true <- log(c(head(th.true, -1), gs.true, s2.true, nu.true - 1))
gamma.true <- tail(th.true, 1)

## FORK CLUSTER
# cl <- makeCluster(nCores, type = "FORK")
# clusterSetRNGStream(cl, 123)  ## set seed
# clusterEvalQ(cl, sink(paste0(printPath, Sys.getpid(), ".txt")))

## PSOCK CLUSTER
cl <- makeCluster(nCores, type = "PSOCK")
clusterSetRNGStream(cl, 123)  ## set seed
clusterExport(cl=cl, as.list(setdiff(ls(), "cl")),
              envir=environment())
clusterEvalQ(cl, sink(paste0(printPath, Sys.getpid(), ".txt")))
clusterEvalQ(cl, library("stringr"))
clusterEvalQ(cl, library("pracma"))
clusterEvalQ(cl, library("Deriv"))
clusterEvalQ(cl, library("deSolve"))
clusterEvalQ(cl, library("Matrix"))
clusterEvalQ(cl, library("mvtnorm"))
clusterEvalQ(cl, library("tmvtnorm"))
clusterEvalQ(cl, library("parallel"))
clusterEvalQ(cl, library("nloptr"))
clusterEvalQ(cl, library("spsUtil"))
clusterEvalQ(cl, library("lbfgsb3c"))
clusterEvalQ(cl, library("optimx"))

Y0s <- parSapply(cl = cl, X = 1:nSim, FUN = function(j){
  return(get.x.em(theta = th.true,
                  nu = nu.true,
                  times = tps,
                  chs.lst = chs))
}, simplify = "array")
stopCluster(cl)

X0s <- get.logistic(Y0s)

# cl <- makeCluster(nCores, type = "FORK")
## FORK CLUSTER

## PSOCK CLUSTER
cl <- makeCluster(nCores, type = "PSOCK")
clusterSetRNGStream(cl, 123)  ## set seed
clusterExport(cl=cl, as.list(setdiff(ls(), "cl")),
              envir=environment())
clusterEvalQ(cl, sink(paste0(printPath, Sys.getpid(), ".txt")))
clusterEvalQ(cl, library("stringr"))
clusterEvalQ(cl, library("pracma"))
clusterEvalQ(cl, library("Deriv"))
clusterEvalQ(cl, library("deSolve"))
clusterEvalQ(cl, library("Matrix"))
clusterEvalQ(cl, library("mvtnorm"))
clusterEvalQ(cl, library("tmvtnorm"))
clusterEvalQ(cl, library("parallel"))
clusterEvalQ(cl, library("nloptr"))
clusterEvalQ(cl, library("spsUtil"))
clusterEvalQ(cl, library("lbfgsb3c"))
clusterEvalQ(cl, library("optimx"))

y0s <- parSapply(cl = cl, X = 1:nSim, FUN = function(j){
  return(get.y.em(times = tps,
                  X = X0s[,,j],
                  gs = gs.true,
                  nu = nu.true,
                  E = E.true,
                  s2 = s2.true))
}, simplify = "array")
stopCluster(cl)
gc()

tps.all <- tps
y0s.all <- y0s
Vsin.all <- Vsin

# sIndex <- sort(sample(1:length(tps.all), length(tps.all)/100*50))

if(pFreq != 100){
  sIndex <- c(1, sort(sample(2:(length(tps.all)-1), round((length(tps.all) - 2)/100*pFreq))), length(tps.all))
  tps <- tps.all[sIndex]
  y0s <- y0s.all[sIndex,]
  Vsin <- Vsin.all[sIndex]
}



cat(paste0("\t Fitting ODEs...\n"))

psi0 <- c(log(rep(1e-2, length(JSENS))),
          log(1),
          log(1e-3))

## FORK CLUSTER
# cl <- makeCluster(nCores, type = "FORK")
# clusterEvalQ(cl, sink(paste0(printPath, Sys.getpid(), ".txt")))

## PSOCK CLUSTER
cl <- makeCluster(nCores, type = "PSOCK")
clusterExport(cl=cl, as.list(setdiff(ls(), "cl")),
              envir=environment())
clusterEvalQ(cl, sink(paste0(printPath, Sys.getpid(), ".txt")))
clusterEvalQ(cl, library("stringr"))
clusterEvalQ(cl, library("pracma"))
clusterEvalQ(cl, library("Deriv"))
clusterEvalQ(cl, library("deSolve"))
clusterEvalQ(cl, library("Matrix"))
clusterEvalQ(cl, library("mvtnorm"))
clusterEvalQ(cl, library("tmvtnorm"))
clusterEvalQ(cl, library("parallel"))
clusterEvalQ(cl, library("nloptr"))
clusterEvalQ(cl, library("spsUtil"))
clusterEvalQ(cl, library("lbfgsb3c"))
clusterEvalQ(cl, library("optimx"))


resODE.all <- parLapply(cl = cl, X = 1:nSim, fun = function(j){
  
  t.proc <- system.time(res.j <- try(optim(par = psi0,
                                           fn = nl_try,
                                           gr = gnl_try,
                                           p9 = tail(th.true,1),
                                           y = y0s[,j],
                                           times = tps,
                                           V = Vsin,
                                           E = E.true,
                                           hmax = 1, # max(tps)
                                           method = "L-BFGS-B",
                                           control = list(lmm = 250,
                                                          factr = 0,
                                                          maxit = 10000,
                                                          trace = 1,
                                                          REPORT = 1)), 
                                     silent = TRUE))
  
  res.j$t.proc <- t.proc
  
  return(res.j)})
stopCluster(cl)
gc()

cat(paste0("\t Fitting nu and s2...\n"))

## FORK CLUSTER
# cl <- makeCluster(nCores, type = "FORK")
# clusterEvalQ(cl, sink(paste0(printPath, Sys.getpid(), ".txt")))

## PSOCK CLUSTER
cl <- makeCluster(nCores, type = "PSOCK")
clusterExport(cl=cl, as.list(setdiff(ls(), "cl")),
              envir=environment())
clusterEvalQ(cl, sink(paste0(printPath, Sys.getpid(), ".txt")))
clusterEvalQ(cl, library("stringr"))
clusterEvalQ(cl, library("pracma"))
clusterEvalQ(cl, library("Deriv"))
clusterEvalQ(cl, library("deSolve"))
clusterEvalQ(cl, library("Matrix"))
clusterEvalQ(cl, library("mvtnorm"))
clusterEvalQ(cl, library("tmvtnorm"))
clusterEvalQ(cl, library("parallel"))
clusterEvalQ(cl, library("nloptr"))
clusterEvalQ(cl, library("spsUtil"))
clusterEvalQ(cl, library("lbfgsb3c"))
clusterEvalQ(cl, library("optimx"))

res.s2nu.all <- parLapply(cl = cl, X = 1:nSim, fun = function(j){
  
  th.j <- exp(head(resODE.all[[j]]$par, 8))
  
  xP.wg.true_j <- try(vode(y = c(get.xP.steady(t = max(tps),
                                              parms = c(th.j, gamma.true)),
                                 get.xP.steady_djs(t = max(tps),
                                                  parms = c(th.j, gamma.true))),
                           times = tps,
                           func = get.dxdP.wg,
                           jacfunc = get.Jf.xP.wg,
                           mf = -21,
                           hmax = 1, # max(tps)
                           parms = c(th.j, gamma.true),
                           rtol = 1e-6),
                      silent = TRUE)[,-1]
  
  psi0j <- c(resODE.all[[j]]$par[10], log(49))
  
  t.proc <- system.time(res.j <- try(optim(par = psi0j, # psi0
                                           fn = nl.nus2.try,
                                           gr = gnl.nus2.try,
                                           g = exp(resODE.all[[j]]$par[9]),
                                           gamma = gamma.true,
                                           times = tps,
                                           V = Vsin,
                                           E = E.true,
                                           y = y0s[,j],
                                           xP.wg = xP.wg.true_j,
                                           method = "L-BFGS-B",
                                           control = list(lmm = 250,
                                                          factr = 0,
                                                          maxit = 10000,
                                                          trace = 1,
                                                          REPORT = 1)), 
                                     silent = TRUE))
  
  res.j$t.proc <- t.proc
  
  return(res.j)})
stopCluster(cl)
gc()

cat(paste0("\t Fitting SDEs...\n"))

## FORK CLUSTER
# cl <- makeCluster(nCores, type = "FORK")
# clusterEvalQ(cl, sink(paste0(printPath, Sys.getpid(), ".txt")))

## PSOCK CLUSTER
cl <- makeCluster(nCores, type = "PSOCK")
clusterExport(cl=cl, as.list(setdiff(ls(), "cl")),
              envir=environment())
clusterEvalQ(cl, sink(paste0(printPath, Sys.getpid(), ".txt")))
clusterEvalQ(cl, library("stringr"))
clusterEvalQ(cl, library("pracma"))
clusterEvalQ(cl, library("Deriv"))
clusterEvalQ(cl, library("deSolve"))
clusterEvalQ(cl, library("Matrix"))
clusterEvalQ(cl, library("mvtnorm"))
clusterEvalQ(cl, library("tmvtnorm"))
clusterEvalQ(cl, library("parallel"))
clusterEvalQ(cl, library("nloptr"))
clusterEvalQ(cl, library("spsUtil"))
clusterEvalQ(cl, library("lbfgsb3c"))
clusterEvalQ(cl, library("optimx"))

res.all <- parLapply(cl = cl, X = 1:nSim, fun = function(j){
  
  psi0_j <- c(head(resODE.all[[j]]$par, 8),
              resODE.all[[j]]$par[9] - log(1 + exp(res.s2nu.all[[j]]$par[2])),
              res.s2nu.all[[j]]$par[1],
              res.s2nu.all[[j]]$par[2])
  
  t.proc <- system.time(res.j <- try(optim(par = psi0_j, # psi0
                                           fn = nl.try,
                                           gr = gnl.try,
                                           gamma = gamma.true,
                                           times = tps,
                                           V = Vsin,
                                           E = E.true,
                                           y = y0s[,j],
                                           hmax = 1, # max(tps)
                                           method = "L-BFGS-B",
                                           control = list(lmm = 250,
                                                          factr = 0,
                                                          maxit = 10000,
                                                          trace = 1,
                                                          REPORT = 1)),
                                     silent = TRUE))
  
    res.j$t.proc <- t.proc
  
  return(res.j)})
stopCluster(cl)
gc()

cat(paste0("\t Fitting SDEs with EKF...\n"))

## FORK CLUSTER
# cl <- makeCluster(nCores, type = "FORK")
# clusterEvalQ(cl, sink(paste0(printPath, Sys.getpid(), ".txt")))

## PSOCK CLUSTER
cl <- makeCluster(nCores, type = "PSOCK")
clusterExport(cl=cl, as.list(setdiff(ls(), "cl")),
              envir=environment())
clusterEvalQ(cl, sink(paste0(printPath, Sys.getpid(), ".txt")))
clusterEvalQ(cl, library("stringr"))
clusterEvalQ(cl, library("pracma"))
clusterEvalQ(cl, library("Deriv"))
clusterEvalQ(cl, library("deSolve"))
clusterEvalQ(cl, library("Matrix"))
clusterEvalQ(cl, library("mvtnorm"))
clusterEvalQ(cl, library("tmvtnorm"))
clusterEvalQ(cl, library("parallel"))
clusterEvalQ(cl, library("nloptr"))
clusterEvalQ(cl, library("spsUtil"))
clusterEvalQ(cl, library("lbfgsb3c"))
clusterEvalQ(cl, library("optimx"))

res.all.ekf <- parLapply(cl = cl, X = 1:nSim, fun = function(j){
  
  psi0_j <- res.all[[j]]$par
  
  t.proc <- system.time(res.j <- try(optim(par = psi0_j, # psi0
                                           fn = nl.ekf.try,
                                           gr = gnl.ekf.try,
                                           gamma = gamma.true,
                                           times = tps,
                                           V = Vsin,
                                           E = E.true,
                                           y = y0s[,j],
                                           hmax = 1, # max(tps)
                                           method = "L-BFGS-B",
                                           control = list(lmm = 250,
                                                          factr = 0,
                                                          maxit = 10000,
                                                          trace = 1,
                                                          REPORT = 1)),
                                     silent = TRUE))
  
  res.j$t.proc <- t.proc
  
  return(res.j)})
stopCluster(cl)
gc()

save.image(paste0(resPath,"/test-ekf-sim-fit.RData"))


write.csv(x = t(sapply(res.all, function(res){
  res$par
})),
file = paste0(resPath, "test-ekf-theta-SDE_", format(s2.true, scientific = FALSE), "-", nu.true, "-pFreq", pFreq, ".csv"))

write.csv(x = t(sapply(resODE.all, function(res){
  res$par
})),
file = paste0(resPath, "test-ekf-theta-ODE_", format(s2.true, scientific = FALSE), "-", nu.true, "-pFreq", pFreq, ".csv"))

write.csv(x = t(sapply(res.all.ekf, function(res){
  res$par
})),
file = paste0(resPath, "test-ekf-theta-SDEwithEKF_", format(s2.true, scientific = FALSE), "-", nu.true, "-pFreq", pFreq, ".csv"))

mus <- get.pred(psi = res.all[[1]]$par, # apply(sapply(res.all, function(res){res$par}), 1, mean),
                gamma = gamma.true,
                times = tps,
                V = Vsin,
                E = E.true,
                hmax = max(tps))

mu_psd <- (mus[1:length(Vsin) %% 10 == 1, "mu"] + qnorm(1 - 0.05/2)*sqrt(mus[1:length(Vsin) %% 10 == 1, "s"]))
mu_msd <- (mus[1:length(Vsin) %% 10 == 1, "mu"] - qnorm(1 - 0.05/2)*sqrt(mus[1:length(Vsin) %% 10 == 1, "s"]))

mus.ekf <- get.pred.ekf(psi = res.all.ekf[[1]]$par, #apply(sapply(res.all.ekf, function(res){res$par}), 1, mean),
                        gamma = gamma.true,
                        times = tps,
                        V = Vsin,
                        E = E.true,
                        hmax = max(tps),
                        y = y0s[,1])

mu_psd.ekf <- (mus.ekf[1:length(Vsin) %% 10 == 1, "mu"] + qnorm(1 - 0.05/2)*sqrt(mus.ekf[1:length(Vsin) %% 10 == 1, "s"]))
mu_msd.ekf <- (mus.ekf[1:length(Vsin) %% 10 == 1, "mu"] - qnorm(1 - 0.05/2)*sqrt(mus.ekf[1:length(Vsin) %% 10 == 1, "s"]))


pdf(paste0(allFigsPath, "fig9a-pfreq", pFreq, ".pdf"), width = 12, height = 4)
par(mar = c(2,5,.5,.5), mfrow = c(1,1))
###
plot(1, 
     type="n", 
     xlab = "time (s)",
     ylab = "current (nA)",
     xlim=c(0, 8), 
     ylim=range(c(mus[1:length(Vsin) %% 10 == 1, "mu"] - qnorm(1 - 0.05/2)*sqrt(mus[1:length(Vsin) %% 10 == 1, "s"]), 
                    mus[1:length(Vsin) %% 10 == 1, "mu"] + qnorm(1 - 0.05/2)*sqrt(mus[1:length(Vsin) %% 10 == 1, "s"]))),
     cex.axis = 2,
     cex.lab = 2)
polygon(c(tps[1:length(tps) %% 10 == 1]/1000, 
          rev(tps[1:length(tps) %% 10 == 1]/1000)), 
        c(mu_psd, rev(mu_msd)),
        col = scales::alpha("blue",.2), lty = 0)

polygon(c(tps[1:length(tps) %% 10 == 1]/1000, 
          rev(tps[1:length(tps) %% 10 == 1]/1000)), 
        c(mu_psd.ekf, rev(mu_msd.ekf)),
        col = scales::alpha("red",.2), lty = 0)

lines(tps[1:length(tps) %% 10 == 1]/1000,
     y0s[1:length(tps) %% 10 == 1, 1], 
     lty = 1,
     lwd = .5,
     cex.axis = 2,
     cex.lab = 2,
     col = scales::alpha("darkgrey",.5))
lines(tps[1:length(tps) %% 10 == 1]/1000, 
      mus.ekf[1:nrow(mus.ekf) %% 10 == 1, 1], 
      col = scales::alpha("red",.5), 
      lwd = .5)
lines(tps[1:length(tps) %% 10 == 1]/1000, 
      mus[1:nrow(mus) %% 10 == 1, 1], 
      col = 'blue', 
      lwd = .5)

legend(x = "bottomright",
       legend = c("simulated current",
                  "DME fit",
                  "EKF fit"),
       col = c("darkgrey",
               "blue",
               "red"),
       lwd = 10,
       horiz = FALSE,
       bty = 'n',
       cex = 1.8)
dev.off()


pdf(paste0(allFigsPath, "fig9b-pfreq", pFreq, ".pdf"), width = 8, height = 4)
par(mar = c(2,5,.5,.5), mfrow = c(1,1))
###
plot(1, 
     type="n", 
     xlab = "time (s)",
     ylab = "current (nA)",
     xlim = c(3.5, 6.3),
     ylim = c(min(mu_msd[3500:6300]), 
              max(mu_psd[3500:6300])),
     cex.axis = 2,
     cex.lab = 2)
polygon(c(tps[1:length(tps) %% 10 == 1]/1000, 
          rev(tps[1:length(tps) %% 10 == 1]/1000)), 
        c(mu_psd, rev(mu_msd)),
        col = scales::alpha("blue",.2), lty = 0)
polygon(c(tps[1:length(tps) %% 10 == 1]/1000, 
          rev(tps[1:length(tps) %% 10 == 1]/1000)), 
        c(mu_psd.ekf, rev(mu_msd.ekf)),
        col = scales::alpha("red",.2), lty = 0)

lines(tps[1:length(tps) %% 10 == 1]/1000,
     y0s[1:length(tps) %% 10 == 1, 1], 
     lty = 1,
     lwd = .5,
     cex.axis = 2,
     cex.lab = 2,
     col = scales::alpha("darkgrey",.5))
lines(tps[1:length(tps) %% 10 == 1]/1000, 
      mus.ekf[1:nrow(mus.ekf) %% 10 == 1, 1], 
      col = scales::alpha("red",.7), 
      lwd = .5)
lines(tps[1:length(tps) %% 10 == 1]/1000, 
      mus[1:nrow(mus) %% 10 == 1, 1], 
      col = 'blue', 
      lwd = .5)
dev.off()

##########################

SDE.times <- rowSums(cbind(sapply(resODE.all, function(resODE){resODE$t.proc[3]}),
                           sapply(res.s2nu.all, function(res.s2nu){res.s2nu$t.proc[3]}),
                           sapply(res.all, function(res.SDE){res.SDE$t.proc[3]})))/60/60


EKF.times <- rowSums(cbind(sapply(resODE.all, function(resODE){resODE$t.proc[3]}),
                           sapply(res.s2nu.all, function(res.s2nu){res.s2nu$t.proc[3]}),
                           sapply(res.all, function(res.SDE){res.SDE$t.proc[3]}),
                           sapply(res.all.ekf, function(res.EKF){res.EKF$t.proc[3]})))/60/60


norm2.ekf <- apply(sapply(res.all.ekf, function(res.EKF){res.EKF$par}), 2, 
                   function(x){norm(x - psi.true, type = "2")/norm(psi.true, type = "2")})

norm2.sde <- apply(sapply(res.all, function(res.SDE){res.SDE$par}), 2, 
                   function(x){norm(x - psi.true, type = "2")/norm(psi.true, type = "2")})


pdf(paste0(allFigsPath, "fig9d-pfreq", pFreq, ".pdf"), width = 5, height = 4)
par(mar = c(3,5.5,.5,.5))
boxplot(cbind(DME = norm2.sde, 
              EKF = norm2.ekf),
        lwd = 2,
        cex.axis = 2,
        cex.lab = 2,
        ylab = expression("\u2016"~hat(psi) - psi[true]~"\u2016"/"\u2016"~psi[true]~"\u2016"))
        # ylab = expression("||"~hat(psi) - psi[true]~"||"/"||"~psi[true]~"||"))
dev.off()

pdf(paste0(allFigsPath, "fig9e-pfreq", pFreq, ".pdf"), width = 5, height = 4)
par(mar = c(3,5.5,.5,.5))
boxplot(cbind(DME = SDE.times, 
              EKF = EKF.times),
        lwd = 2,
        cex.axis = 2,
        cex.lab = 2,
        ylab = "computing time (hours)")
dev.off()

round(mean(SDE.times),2)
round(mean(EKF.times))



##########################


# pdf(paste0(figPath, "test-ekf-UQ-pFreq", pFreq, ".pdf"), width = 8, height = 4)
# par(mar = c(3,5,.5,.5), mfrow = c(1,1))
# boxplot(t(sapply(res.all, function(res){res$par})/psi.true),
#         cex.axis = 1.5,
#         cex.lab = 1.5,
#         pch = 20,
#         lwd = 2,
#         outline = FALSE,
#         xaxt = 'n')
# abline(h = 1, lwd = 3, col = "grey")
# axis(side = 1,
#      at = 1:length(psi.true),
#      labels = head(c(sapply(1:8, function(j){bquote(theta[.(j)])}), 
#                      expression(g[s]),
#                      expression(sigma^2),
#                      expression(eta),
#                      expression(1)), -1),
#      padj = .5,
#      cex.axis = 2)
# dev.off()


pdf(paste0(allFigsPath, "fig9c-pfreq", pFreq, ".pdf"), width = 6, height = 5)
par(mar = c(4,5.5,.5,.5))
plot(rep(psi.true, times = nSim), c(sapply(res.all, function(res){
  res$par
})),
xlim = range(rep(psi.true, times = nSim), c(sapply(res.all, function(res){
  res$par
}))),
ylim = range(rep(psi.true, times = nSim), c(sapply(res.all, function(res){
  res$par
}))),
pch = 21,
cex = 2,
cex.axis = 2,
cex.lab = 2,
xlab = expression(psi[true]),
ylab = expression(hat(psi)), 
col = "blue",
lwd = 2)
points(rep(psi.true, times = nSim),
       c(sapply(res.all.ekf, function(res){
  res$par
})),
pch = 21,
cex = 2,
col = "red",
lwd = 2)
lines(range(rep(psi.true, times = nSim), c(sapply(res.all, function(res){
  res$par
}))),
range(rep(psi.true, times = nSim), c(sapply(res.all, function(res){
  res$par
}))),
lwd = 5,
col = "grey")
legend(x = "bottomright",
       legend = c("DME",
                  "EKF"),
       col = c("blue", "red"),
       bty ="n",
       cex = 2,
       lwd = 10)
dev.off()

