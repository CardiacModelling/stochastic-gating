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
            "nloptr")
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
source("./source/sim-fit-nCh-fun.R")

# s2_nu.true.min.max_nSim_nCores <- commandArgs(trailingOnly = TRUE) # s2_nu.true.min.max_nSim_nCores <- "0.00001-100-10000-10-10"
s2.true <- as.numeric(unlist(strsplit(s2_nu.true.min.max_nSim_nCores, split = "-"))[1])
nu.true.min <- as.numeric(unlist(strsplit(s2_nu.true.min.max_nSim_nCores, split = "-"))[2])
nu.true.max <- as.numeric(unlist(strsplit(s2_nu.true.min.max_nSim_nCores, split = "-"))[3])
nSim <- as.numeric(unlist(strsplit(s2_nu.true.min.max_nSim_nCores, split = "-"))[4])
nCores <- as.numeric(unlist(strsplit(s2_nu.true.min.max_nSim_nCores, split = "-"))[5])

cat(paste0("\ts2 = ", s2.true, 
           "\tnu.min = ", nu.true.min,
           "\tnu.max = ", nu.true.max,
           "\tnSim = ", nSim, 
           "\tnCores = ", nCores, "\n"))

##############
## currTime ##
##############

# currTime <- print(as.numeric(Sys.time())*1000, digits=20)
currTime <- Sys.time()
currTime <- gsub(" ", "_", currTime)
currTime <- gsub("-", "", currTime)
currTime <- gsub(":", "", currTime)
currTime <- paste(currTime, s2_nu.true.min.max_nSim_nCores, sep = "_")

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
gamma.true <- tail(th.true, 1)


nu.trues <- seq(from = nu.true.min,
                to = nu.true.max,
                length.out = nSim)

cl <- makeCluster(nCores, type = "FORK")
clusterEvalQ(cl, sink(paste0(printPath, Sys.getpid(), ".txt")))

Y0s <- parSapply(cl = cl, X = 1:nSim, FUN = function(j){
  set.seed(123)
  return(get.x.em(theta = th.true,
                  nu = nu.trues[j],
                  times = tps,
                  chs.lst = chs))
}, simplify = "array")
stopCluster(cl)

X0s <- get.logistic(Y0s)

cl <- makeCluster(nCores, type = "FORK")
y0s <- parSapply(cl = cl, X = 1:nSim, FUN = function(j){
  return(get.y.em(times = tps,
                  X = X0s[,,j],
                  gs = g.true/nu.trues[j],
                  nu = nu.trues[j],
                  E = E.true,
                  s2 = s2.true))
}, simplify = "array")
stopCluster(cl)
gc()

cat(paste0("\t Fitting ODEs...\n"))

psi0 <- c(log(rep(1e-2, length(JSENS))),
          log(1),
          log(1e-3))

cl <- makeCluster(nCores, type = "FORK")
clusterEvalQ(cl, sink(paste0(printPath, Sys.getpid(), ".txt")))

resODE.all <- parLapply(cl = cl, X = 1:nSim, fun = function(j){
  
  res.j <- try(optim(par = psi0,
                     fn = nl_try,
                     gr = gnl_try,
                     p9 = tail(th.true,1),
                     y = y0s[,j],
                     times = tps,
                     V = Vsin,
                     E = E.true,
                     hmax = 1,
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

cat(paste0("\t Fitting nu and s2...\n"))

cl <- makeCluster(nCores, type = "FORK")
clusterEvalQ(cl, sink(paste0(printPath, Sys.getpid(), ".txt")))

res.s2nu.all <- parLapply(cl = cl, X = 1:nSim, fun = function(j){
  
  th.j <- exp(head(resODE.all[[j]]$par, 8))
  
  xP.wg.true_j <- try(vode(y = c(get.x.steady(t = max(tps),
                                              parms = c(th.j, gamma.true)),
                                 rep(0, 10), # nStates*(nStates - 1)/2 + nStates
                                 get.x.steady_djs(t = max(tps),
                                                  parms = c(th.j, gamma.true)),
                                 rep(0, 80)),
                           times = tps,
                           func = get.dxdP.wg,
                           jacfunc = get.Jf.xP.wg,
                           mf = -21,
                           hmax = 1, # max(tps),
                           parms = c(th.j, gamma.true),
                           rtol = 1e-6),
                      silent = TRUE)[,-1]
  
  psi0j <- c(resODE.all[[j]]$par[10], log(49))
  
  res.j <- try(optim(par = psi0j, # psi0
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
               silent = TRUE)
  
  return(res.j)})
stopCluster(cl)
gc()

cat(paste0("\t Fitting SDEs...\n"))

cl <- makeCluster(nCores, type = "FORK")
clusterEvalQ(cl, sink(paste0(printPath, Sys.getpid(), ".txt")))

res.all <- parLapply(cl = cl, X = 1:nSim, fun = function(j){
  
  psi0_j <- c(head(resODE.all[[j]]$par, 8),
              resODE.all[[j]]$par[9] - log(1 + exp(res.s2nu.all[[j]]$par[2])),
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
                     hmax = 1, # max(tps),
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

cat(paste0("\t Estimating n. of channels with ODEs...\n"))
cl <- makeCluster(nCores, type = "FORK")
clusterEvalQ(cl, sink(paste0(printPath, Sys.getpid(), ".txt")))

res.nchODEs.all <- parLapply(cl = cl, X = 1:nSim, fun = function(j){
  
  theta.pred <- c(exp(head(resODE.all[[j]]$par, -2)),
                  tail(th.true, 1))
  
  system.time(x.pred <- try(vode(y = c(get.x.steady(t = max(tps),
                                                    parms = theta.pred),
                                       get.x.steady_djs(t = max(tps),
                                                        parms = theta.pred)),
                                 times = tps,
                                 func = get.dx.wg,
                                 jacfunc = get.Jf.wg,
                                 mf = -21,
                                 hmax = 1,
                                 parms = theta.pred), silent = TRUE))
  
  mu.ODE <- x.pred[,"O"]*exp(tail(resODE.all[[j]]$par, 2)[1])*(Vsin - E.true)
  s2.ODE <- rep(exp(tail(resODE.all[[j]]$par, 1)), length(mu.ODE)) ## iid noise estimated by ODEs
  
  n2c.resj <- optim(par = c(1,1),
                  fn = n2c.nl.comp,
                  gr = n2c.gnl.comp,
                  mu = mu.ODE,
                  s2 = s2.ODE,
                  V = Vsin,
                  E = E.true,
                  method = "L-BFGS-B",
                  control = list(lmm = 250,
                                 factr = 0,
                                 trace = 1,
                                 pgtol = 0,
                                 maxit = 10000))
  
  return(n2c.resj)})
stopCluster(cl)
gc()

save.image(paste0(resPath,"/sim-fit-n2c.RData"))

write.csv(x = t(sapply(res.all, function(res){
  res$par
})),
file = paste0(resPath, "sim-fit-theta-SDE.csv"))

write.csv(x = t(sapply(resODE.all, function(res){
  res$par
})),
file = paste0(resPath, "sim-fit-theta-ODE.csv"))


h1 <- hist(sapply(resODE.all, FUN = function(res){
  tail(res$par, 1)
}), plot = FALSE)

h2 <- hist(sapply(res.all, FUN = function(res){
  tail(res$par, 2)[1]
}), plot = FALSE)

pdf(paste0(allFigsPath, "fig6b",as.character(s2.true),".pdf"), width = 8, height = 4)
par(mar = c(3,3,1,1), mfrow = c(1,2))
plot(h1, 
     col = scales::alpha('blue', alpha = .2), 
     lty = 0, 
     xlab = '', 
     freq = TRUE, 
     main = "", 
     cex.axis = 1.2, 
     cex.lab = 1.2,
     ylim = c(0,max(h1$counts, h2$counts)),
     ylab = '')  # first histogram
abline(v = log(s2.true), lty = 2, lwd = 3)
par(mar = c(3,1,1,1))
plot(h2, 
     col = scales::alpha('red', alpha = .2), 
     lty = 0, 
     xlab = '', 
     freq = TRUE, 
     main = "", 
     cex.axis = 1.2, 
     cex.lab = 1.2,
     ylim = c(0, max(h1$counts, h2$counts)),
     yaxt = 'n',
     ylab = '')  # second
abline(v = log(s2.true), lty = 2, lwd = 3)
dev.off()


pdf(paste0(allFigsPath, "fig6a", as.character(s2.true),".pdf"), width = 5, height = 4.5)
par(mar = c(4,5.5,1,1))
plot(log(nu.trues),
     sapply(res.all, function(res){
       tail(res$par,2)
     })[2,],
     ylim = range(c(sapply(res.all, function(res){
       tail(res$par,2)
     })[2,], 
     sapply(res.nchODEs.all, function(res){
       res$par
     })[2,])),
     xlab = expression(log~eta),
     ylab = expression(log~hat(eta)),
     cex.axis = 2,
     cex.lab = 2,
     cex = 2,
     lwd = 2,
     pch = 20)
lines(range(c(log(nu.trues),
     sapply(res.all, function(res){
  tail(res$par,2)
})[2,])),
range(c(log(nu.trues),
     sapply(res.all, function(res){
  tail(res$par,2)
})[2,])), lwd = 5, col = scales::alpha("red", .9))

points(log(nu.trues),
     sapply(res.nchODEs.all, function(res){
  res$par
})[2,], 
pch = 20, 
cex = 2,
col = "blue")
legend(x = 'bottomright',
       legend = c("Gray et al.",
                  "our method",
                  "true"),
       lty = 1,
       lwd = c(NA, NA, 5),
       pch = c(20, 20, NA),
       col = c("blue", "black", "red"),
       cex = 1.5,
       bty = 'n',
       pt.cex = 3,
       pt.lwd = 3,
       text.font = 2)
dev.off()
