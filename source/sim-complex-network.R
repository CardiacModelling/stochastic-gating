
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
source("./source/sim-complex-network-fun.R")

# s2_nu.true_nSim_nCores <- commandArgs(trailingOnly = TRUE) # s2_nu.true_nSim_nCores <- "0.00001-10000-100-30"
s2.true <- as.numeric(unlist(strsplit(s2_nu.true_nSim_nCores, split = "-"))[1])
nu.true <- as.numeric(unlist(strsplit(s2_nu.true_nSim_nCores, split = "-"))[2])
nSim <- as.numeric(unlist(strsplit(s2_nu.true_nSim_nCores, split = "-"))[3])
nCores <- as.numeric(unlist(strsplit(s2_nu.true_nSim_nCores, split = "-"))[4])

cat(paste0("\ts2 = ", s2.true, 
           "\tnu = ", nu.true, 
           "\tnSim = ", nSim, 
           "\tnCores = ", nCores, "\n"))

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

## 8 STATE MODEL
chs <- c("CCC", "ICIC", "CC", "C", "O", "F", "I", "IC")

rcs <- c("C->O",
         "O->C",
         "C->F",
         "F->C",
         "O->I",
         "I->O",
         "F->I",
         "I->F",
         "CC->C",
         "C->CC",
         "IC->I",
         "I->IC",
         "C->IC",
         "IC->C",
         "ICIC->IC",
         "IC->ICIC",
         "CC->ICIC",
         "ICIC->CC",
         "CCC->CC",
         "CC->CCC")

cnstrs <- c("lambda\\[\\'IC->I\\'\\]=lambda\\[\\'C->O\\'\\]",
            "lambda\\[\\'I->IC\\'\\]=lambda\\[\\'O->C\\'\\]",
            "lambda\\[\\'C->IC\\'\\]=lambda\\[\\'O->I\\'\\]",
            "lambda\\[\\'IC->C\\'\\]=lambda\\[\\'I->O\\'\\]",
            "lambda\\[\\'F->C\\'\\]=lambda\\[\\'O->C\\'\\]",
            "lambda\\[\\'F->I\\'\\]=lambda\\[\\'O->I\\'\\]",
            "lambda\\[\\'C->F\\'\\]=lambda\\[\\'C->O\\'\\]",
            "lambda\\[\\'I->F\\'\\]=lambda\\[\\'I->O\\'\\]",
            ##
            "lambda\\[\\'ICIC->IC\\'\\]=lambda\\[\\'CC->C\\'\\]",
            "lambda\\[\\'IC->ICIC\\'\\]=lambda\\[\\'C->CC\\'\\]",
            "lambda\\[\\'CC->ICIC\\'\\]=lambda\\[\\'O->I\\'\\]",
            "lambda\\[\\'ICIC->CC\\'\\]=lambda\\[\\'I->O\\'\\]")


# ## synthetic trace from 5-state model
psi.true.ODE <- c(-2.2002062770407952996265521505847573280334472656250000,
                  -6.3831902530543862539502697472926229238510131835937500,
                  -5.3795985179253715813274538959376513957977294921875000,
                  -3.9087510966015930335970551823265850543975830078125000,
                  -2.2655636546359740890466127893887460231781005859375000,
                  -4.6357680632051785352132355910725891590118408203125000,
                  -5.2448566783358661069769368623383343219757080078125000,
                  -3.4661739234114912200368507910752668976783752441406250,
                  -3.6440308727120389598042038414860144257545471191406250,
                  -6.1738714124377853664782378473319113254547119140625000,
                  -7.4630360734451075543915976595599204301834106445312500,
                  -3.1373606866316157137930531462188810110092163085937500,
                  -8.0783750124985953533496285672299563884735107421875000,
                  -2.7254649246451005950575563474558293819427490234375000,
                  -3.6387576283646549946126924623968079686164855957031250,
                  -3.3370094437059600878114906663540750741958618164062500,
                  -1.7532388960943432465455771307460963726043701171875000,
                  -9.8977136306510420382664960925467312335968017578125000)

th.true <- c(exp(head(psi.true.ODE, 16)), c(1/6.7)/sum(c(1/6.7,
                                                         1/2.5)))

g.true <- exp(psi.true.ODE[17])

FLKR <- TRUE
ntheta <- (length(rcs) - length(cnstrs))*2 + 1
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
# E.true <- -88.35746
E.true <- -88.3574598824808532526731141842901706695556640625
gs.true <- g.true/nu.true

psi.true <- log(c(head(th.true, -1), gs.true, s2.true, nu.true - 1))
gamma.true <- tail(th.true, 1)

# # FORK CLUSTER
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

# ## FORK CLUSTER
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


cat(paste0("\t Fitting ODEs...\n"))

psi0 <- c(log(rep(1e-2, length(JSENS))),
          log(1),
          log(1e-3))

# # FORK CLUSTER
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
  
  res.j <- try(optim(par = psi0,
                     fn = nl_try,
                     gr = gnl_try,
                     p17 = tail(th.true,1),
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
               silent = TRUE)
  
  # res.j <- try(parma::cmaes(pars = res.j$par, # resODE$par
  #                           fun = nl_try,
  #                           lower = rep(-Inf, length(psi0)),
  #                           upper = rep(Inf, length(psi0)),
  #                           insigma = .1,
  #                           ctrl = CMAES.control,
  #                           p17 = gamma.true,
  #                           y = y0s[,j],
  #                           times = tps,
  #                           V = Vsin,
  #                           E = E.true,
  #                           hmax = 1),
  #              silent = TRUE)
  
  return(res.j)})
stopCluster(cl)
gc()

cat(paste0("\t Fitting nu and s2...\n"))

# # FORK CLUSTER
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
  
  th.j <- exp(head(resODE.all[[j]]$par, 16))
  # th.j <- exp(head(resODE.all[[j]]$bestever$x, 16))
  
  # xP.wg.true_j <- try(vode(y = c(get.xP.steady(t = max(tps),
  #                                             parms = c(th.j, gamma.true)),
  #                                get.xP.steady_djs(t = max(tps),
  #                                                 parms = c(th.j, gamma.true))),
  #                          times = tps,
  #                          func = get.dxdP.wg,
  #                          jacfunc = get.Jf.xP.wg,
  #                          mf = -21,
  #                          hmax = 1, # max(tps)
  #                          parms = c(th.j, gamma.true),
  #                          rtol = 1e-6),
  #                     silent = TRUE)[,-1]
  
  xP.wg.true_j <- try(vode(y = c(get.xP.steady(t = max(tps),
                                               parms = c(th.j, gamma.true))),
                           times = tps,
                           func = get.dxdP,
                           jacfunc = get.Jf.xP,
                           mf = -21,
                           hmax = 1, # max(tps)
                           parms = c(th.j, gamma.true),
                           rtol = 1e-6),
                      silent = TRUE)[,-1]
  
  
  psi0j <- c(resODE.all[[j]]$par[18], log(49))
  
  res.j <- try(optim(par = psi0j, # psi0
                     fn = nl.nus2.try,
                     gr = gnl.nus2.try,
                     g = exp(resODE.all[[j]]$par[17]),
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
               silent = TRUE)
  
  return(res.j)})
stopCluster(cl)
gc()

cat(paste0("\t Fitting SDEs...\n"))

# # FORK CLUSTER
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
  
  psi0_j <- c(head(resODE.all[[j]]$par, 16),
              resODE.all[[j]]$par[17] - log(1 + exp(res.s2nu.all[[j]]$par[2])),
              res.s2nu.all[[j]]$par[1],
              res.s2nu.all[[j]]$par[2])
  
  res.j <- try(optim(par = psi0_j, # psi0
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
               silent = TRUE)
  
  return(res.j)})
stopCluster(cl)
gc()

save.image(paste0(resPath,"/complexNet-opt-sim-fit.RData"))


write.csv(x = t(sapply(res.all, function(res){
  res$par
})),
file = paste0(resPath, "complexNet-fit-theta-SDE_", format(s2.true, scientific = FALSE), "-", nu.true, ".csv"))

write.csv(x = t(sapply(resODE.all, function(res){
  res$par
})),
file = paste0(resPath, "complexNet-sim-fit-theta-ODE_", format(s2.true, scientific = FALSE), "-", nu.true, ".csv"))

mus <- get.pred(psi = res.all[[1]]$par,
                gamma = gamma.true,
                times = tps,
                V = Vsin,
                E = E.true,
                hmax = max(tps))
# mus[which(mus[,2] < 0),2] <- exp(tail(res.all[[1]]$par,2)[1])

mu_psd <- (mus[1:length(Vsin) %% 10 == 1, "mu"] + qnorm(1 - 0.05/2)*sqrt(mus[1:length(Vsin) %% 10 == 1, "s"]))
mu_msd <- (mus[1:length(Vsin) %% 10 == 1, "mu"] - qnorm(1 - 0.05/2)*sqrt(mus[1:length(Vsin) %% 10 == 1, "s"]))


pdf(paste0(allFigsPath, "fig8c.pdf"), width = 12, height = 4)
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
###

lines(tps[1:length(tps) %% 10 == 1]/1000,
     y0s[1:length(y0s[,1]) %% 10 == 1,1], 
     lty = 1,
     lwd = .5,
     cex.axis = 2,
     cex.lab = 2,
     col = "darkgrey")
lines(tps[1:length(tps) %% 10 == 1]/1000, 
      mus[1:nrow(mus) %% 10 == 1, 1], 
      col = 'blue', 
      lwd = .5)
legend(x = "bottomright",
       legend = c("simulated current",
                  "first two moments"),
       col = c("darkgrey",
               "blue"),
       lwd = 10,
       horiz = FALSE,
       bty = 'n',
       cex = 1.8)
dev.off()


pdf(paste0(allFigsPath, "fig8c-zoom.pdf"), width = 8, height = 4)
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
###

lines(tps[1:length(tps) %% 10 == 1]/1000,
     y0s[1:nrow(y0s) %% 10 == 1,1], 
     lty = 1,
     lwd = .5,
     cex.axis = 2,
     cex.lab = 2,
     col = "darkgrey")
lines(tps[1:length(tps) %% 10 == 1]/1000, 
      mus[1:nrow(mus) %% 10 == 1, 1], 
      col = 'blue', 
      lwd = .5)
# legend(x = "bottomright",
#        legend = c("simulated current",
#                   "first two moments"),
#        col = c("darkgrey",
#                "blue"),
#        lwd = 10,
#        horiz = FALSE,
#        bty = 'n',
#        cex = 1.8)
dev.off()



# pdf(paste0(figPath, "complexNet-UQ.pdf"), width = 8, height = 4)
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
#      labels = head(c(sapply(1:16, function(j){bquote(theta[.(j)])}), 
#                      expression(g[s]),
#                      expression(sigma^2),
#                      expression(eta),
#                      expression(1)), -1),
#      padj = .5,
#      cex.axis = 2)
# dev.off()


pdf(paste0(allFigsPath, "fig8b.pdf"), width = 6, height = 5)
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
ylab = expression(hat(psi)))
lines(range(rep(psi.true, times = nSim), c(sapply(res.all, function(res){
  res$par
}))),
range(rep(psi.true, times = nSim), c(sapply(res.all, function(res){
  res$par
}))),
lwd = 5,
col = "grey")
dev.off()
