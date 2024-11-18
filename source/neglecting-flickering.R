
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
source("./source/neglecting-flickering-fun.R")

# s2_nu.true <- commandArgs(trailingOnly = TRUE) # s2_nu.true <- "0.00001-1000"
s2.true <- as.numeric(unlist(strsplit(s2_nu.true, split = "-"))[1])
nu.true <- as.numeric(unlist(strsplit(s2_nu.true, split = "-"))[2])

cat(paste0("\ts2 = ", s2.true, 
           "\tnu = ", nu.true, "\n"))

##############
## currTime ##
##############

currTime <- Sys.time()
currTime <- gsub(" ", "_", currTime)
currTime <- gsub("-", "", currTime)
currTime <- gsub(":", "", currTime)
currTime <- paste(currTime, s2_nu.true, sep = "_")

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

######################################
## MODEL WITH FLICKERING MECHANISM
######################################

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


set.seed(123)
Y0 <- get.x.em(theta = th.true,
               nu = nu.true,
               times = tps,
               chs.lst = chs)

X0 <- get.logistic(Y0)
y0 <- get.y.em(times = tps,
               X = X0,
               gs = gs.true,
               nu = nu.true,
               E = E.true,
               s2 = s2.true)

cat(paste0("\t Fitting ODEs for MLE...\n"))

psi0 <- c(log(rep(1e-2, length(JSENS))),
          log(1),
          log(1e-3))

resODE <- try(optim(par = psi0,
                    fn = nl_try,
                    gr = gnl_try,
                    p9 = tail(th.true,1),
                    y = y0,
                    times = tps,
                    V = Vsin,
                    E = E.true,
                    hmax = 1,
                    method = "L-BFGS-B",
                    control = list(lmm = 2500,
                                   factr = 0,
                                   maxit = 10000,
                                   trace = 1,
                                   REPORT = 1)), 
              silent = TRUE)

cat(paste0("\t Fitting nu and s2 for MLE...\n"))

th.mle <- exp(head(resODE$par, 8))

xP.wg.mle <- try(vode(y = c(get.x.steady(t = max(tps),
                                            parms = c(th.mle, gamma.true)),
                               rep(0, 10), # nStates*(nStates - 1)/2 + nStates
                               get.x.steady_djs(t = max(tps),
                                                parms = c(th.mle, gamma.true)),
                               rep(0, 80)),
                         times = tps,
                         func = get.dxdP.wg,
                         jacfunc = get.Jf.xP.wg,
                         mf = -21,
                         hmax = 1, # max(tps),
                         parms = c(th.mle, gamma.true),
                         rtol = 1e-6),
                    silent = TRUE)[,-1]

psi0.mle <- c(resODE$par[10], log(49))

res.s2nu <- try(optim(par = psi0.mle, # psi0
                   fn = nl.nus2.try,
                   gr = gnl.nus2.try,
                   g = exp(resODE$par[9]),
                   gamma = gamma.true,
                   times = tps,
                   V = Vsin,
                   E = E.true,
                   y = y0,
                   xP.wg = xP.wg.mle,
                   method = "L-BFGS-B",
                   control = list(lmm = 250,
                                  factr = 0,
                                  maxit = 10000,
                                  trace = 1,
                                  REPORT = 1)), 
             silent = TRUE)

# exp(res.s2nu$par)

cat(paste0("\t Fitting SDEs for MLE...\n"))
psi0 <- c(head(resODE$par, 8),
            resODE$par[9] - log(1 + exp(res.s2nu$par[2])),
            res.s2nu$par[1],
            res.s2nu$par[2])

res.SDE <- try(optim(par = psi0, # psi0
                     fn = nl.try,
                     gr = gnl.try,
                     gamma = gamma.true,
                     times = tps,
                     V = Vsin,
                     E = E.true,
                     y = y0,
                     hmax = 1, # max(tps),
                     method = "L-BFGS-B",
                     control = list(lmm = 2500,
                                    factr = 0,
                                    maxit = 10000,
                                    trace = 1,
                                    REPORT = 1),
                     hessian = FALSE),
               silent = TRUE)

mus <- get.pred(psi = res.SDE$par,
                gamma = gamma.true,
                times = tps,
                V = Vsin,
                E = E.true,
                hmax = max(tps))

######################################
## MODEL WITHOUT FLICKERING MECHANISM
######################################

chs <- c("C", "O", "I", "IC")
rcs <- c("C->O",
         "O->C",
         "O->I",
         "I->O",
         "IC->I",
         "I->IC",
         "C->IC",
         "IC->C")

cnstrs <- c("lambda\\[\\'IC->I\\'\\]=lambda\\[\\'C->O\\'\\]",
            "lambda\\[\\'I->IC\\'\\]=lambda\\[\\'O->C\\'\\]",
            "lambda\\[\\'C->IC\\'\\]=lambda\\[\\'O->I\\'\\]",
            "lambda\\[\\'IC->C\\'\\]=lambda\\[\\'I->O\\'\\]")

FLKR <- FALSE
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
                  flickering = FLKR)

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


cat(paste0("\t Fitting ODEs for MLE, no flickering mechanism...\n"))

psi0 <- c(log(rep(1e-2, length(JSENS))),
          log(1),
          log(1e-3))

resODE_noF <- try(optim(par = psi0,
                    fn = nl_try,
                    gr = gnl_try,
                    p9 = tail(th.true,1),
                    y = y0,
                    times = tps,
                    V = Vsin,
                    E = E.true,
                    hmax = 1,
                    method = "L-BFGS-B",
                    control = list(lmm = 2500,
                                   factr = 0,
                                   maxit = 10000,
                                   trace = 1,
                                   REPORT = 1)), 
              silent = TRUE)

cat(paste0("\t Fitting nu and s2 for MLE...\n"))

th.mle <- exp(head(resODE_noF$par, 8))

xP.wg.mle_noF <- try(vode(y = c(get.x.steady(t = max(tps),
                                            parms = c(th.mle, gamma.true)),
                               rep(0, 6), # nStates*(nStates - 1)/2 + nStates
                               get.x.steady_djs(t = max(tps),
                                                parms = c(th.mle, gamma.true)),
                               rep(0, 6*8)),
                         times = tps,
                         func = get.dxdP.wg,
                         jacfunc = get.Jf.xP.wg,
                         mf = -21,
                         hmax = 1, # max(tps),
                         parms = c(th.mle, gamma.true),
                         rtol = 1e-6),
                    silent = TRUE)[,-1]

psi0.mle <- c(resODE_noF$par[10], log(49))

res.s2nu_noF <- try(optim(par = psi0.mle, # psi0
                   fn = nl.nus2.try_noF,
                   gr = gnl.nus2.try_noF,
                   g = exp(resODE_noF$par[9]),
                   gamma = gamma.true,
                   times = tps,
                   V = Vsin,
                   E = E.true,
                   y = y0,
                   xP.wg = xP.wg.mle_noF,
                   method = "L-BFGS-B",
                   control = list(lmm = 250,
                                  factr = 0,
                                  maxit = 10000,
                                  trace = 1,
                                  REPORT = 1)), 
             silent = TRUE)



cat(paste0("\t Fitting SDEs for MLE...\n"))
psi0 <- c(head(resODE_noF$par, 8),
            resODE_noF$par[9] - log(1 + exp(res.s2nu_noF$par[2])),
            res.s2nu_noF$par[1],
            res.s2nu_noF$par[2])

res.SDE_noF <- try(optim(par = psi0, # psi0
                     fn = nl.try_noF,
                     gr = gnl.try_noF,
                     gamma = gamma.true,
                     times = tps,
                     V = Vsin,
                     E = E.true,
                     y = y0,
                     hmax = 1, # max(tps),
                     method = "L-BFGS-B",
                     control = list(lmm = 2500,
                                    factr = 0,
                                    maxit = 10000,
                                    trace = 1,
                                    REPORT = 1),
                     hessian = FALSE),
               silent = TRUE)

mus_noF <- get.pred_noF(psi = res.SDE_noF$par,
                        gamma = gamma.true,
                        times = tps,
                        V = Vsin,
                        E = E.true,
                        hmax = max(tps))

(1 + exp(tail(res.SDE_noF$par,1)))/(1 + exp(tail(psi.true,1)))
(1 + exp(tail(res.SDE$par,1)))/(1 + exp(tail(psi.true,1)))
c(1/6.7)/sum(c(1/6.7,
               1/2.5)) # \pi_1
exp(psi.true[9])/exp(res.SDE_noF$par[9])

pdf(paste0(allFigsPath, "fig10a.pdf"), width = 4.5, height = 4.5)
par(mar = c(4,5.5,.5,.5))
plot(x = psi.true, 
     y = res.SDE_noF$par,
     ylim = range(psi.true,
                  res.SDE_noF$par),
     xlim = range(psi.true,
                  res.SDE_noF$par),
     pch = "*",
     cex = 3,
     cex.axis = 2,
     cex.lab = 2,
     col = c(rep("black", 10), "red"),
     xlab = expression(log~psi),
     ylab = expression(log~hat(psi)))
lines(range(psi.true,
                  res.SDE_noF$par),
      range(psi.true,
                  res.SDE_noF$par),
      col = 'grey',
      lwd = 5)
points(x = psi.true, 
       y = res.SDE$par,
            col = c(rep("black", 10), "red"),
       pch = 21,
       lwd = 3,
       cex = 3)
legend(x = "bottomright",
       legend = c("5-state model",
                  "4-state model"),
       lty = NA,
       lwd = 2,
       pch = c(21, 8),
       bty = 'n',
       cex = 1.5)
dev.off()


mu_psd <- (mus[1:length(Vsin) %% 10 == 1, "mu"] + qnorm(1-0.05/2)*sqrt(mus[1:length(Vsin) %% 10 == 1, "s"]))
mu_msd <- (mus[1:length(Vsin) %% 10 == 1, "mu"] - qnorm(1-0.05/2)*sqrt(mus[1:length(Vsin) %% 10 == 1, "s"]))

mu_psd_noF <- (mus_noF[1:length(Vsin) %% 10 == 1, "mu"] + qnorm(1-0.05/2)*sqrt(mus[1:length(Vsin) %% 10 == 1, "s"]))
mu_msd_noF <- (mus_noF[1:length(Vsin) %% 10 == 1, "mu"] - qnorm(1-0.05/2)*sqrt(mus[1:length(Vsin) %% 10 == 1, "s"]))

pdf(paste0(allFigsPath, "fig10b.pdf"), width = 12, height = 4.5)
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
        c(mu_psd_noF, rev(mu_msd_noF)),
        col = scales::alpha("red",.2), lty = 0)
###
lines(x = tps[1:length(tps) %% 10 == 1]/1000, 
      y = y0[1:length(tps) %% 10 == 1], 
      type = 'l', 
      lty = 1,
      lwd = .5,
      col = scales::alpha("black",.5),
      xlab = 'time (s)',
      ylab = 'current',
      cex.axis = 2,
      cex.lab = 2)
###
lines(tps[1:length(tps) %% 10 == 1]/1000, 
      mus[1:length(Vsin) %% 10 == 1, "mu"], 
      col = scales::alpha("blue",.5),
      lwd = .5)
lines(tps[1:length(tps) %% 10 == 1]/1000, 
      mus_noF[1:length(Vsin) %% 10 == 1, "mu"], 
      col = scales::alpha("red",.5),
      lwd = .5)
legend(x = "bottomright",
       legend = c("simulated data",
                  "5-state model",
                  "4-state model"),
       lty = c(1,1),
       lwd = 7,
       col = c("black", "blue", "red"),
       bty = 'n',
       cex = 1.5)
dev.off()
