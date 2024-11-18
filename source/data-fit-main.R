
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
source("./source/data-fit-fun.R")

# nCell <- as.numeric(commandArgs(trailingOnly = TRUE)) # nCell <- "1"

cat(nCell)

##############
## currTime ##
##############

currTime <- Sys.time()
currTime <- gsub(" ", "_", currTime)
currTime <- gsub("-", "", currTime)
currTime <- gsub(":", "", currTime)
currTime <- paste(currTime, nCell, sep = "_")

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
resData <- paste0(resPath, "/FitData/")
ifelse(!dir.exists(file.path(resData)), dir.create(file.path(resData)), FALSE)
resPath <- paste0(paste(resPath, currTime, sep = "/"), "/")
figPath <- paste0(resPath, "figures/")
printPath <- paste0(resPath, "print/")

ifelse(!dir.exists(file.path(resPath)), dir.create(file.path(resPath)), FALSE)
ifelse(!dir.exists(file.path(figPath)), dir.create(file.path(figPath)), FALSE)
ifelse(!dir.exists(file.path(printPath)), dir.create(file.path(printPath)), FALSE)


ty <- read.csv(file = paste0("./data/sine-wave/cell-", nCell, ".csv"))
ty[,"time"] <- round(ty[,"time"], 2)
to_remove <- which(ty[,"time"] %in% c(sapply(c(250,
                                               300,
                                               500,
                                               1500,
                                               2000,
                                               3000,
                                               6500,
                                               7000), function(t){
                                                 t + seq(0,5,by = .1)
                                               })))

ty <- ty[-to_remove,]
tps <- ty[, "time"]

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
E.true <- reversal_potential(temperature(nCell))


gamma.true <- c(1/6.7)/sum(c(1/6.7,
                             1/2.5))


cat(paste0("\t Fitting ODEs for MLE...\n"))

psi0 <- c(log(rep(1e-2, length(JSENS))),
          log(1),
          log(1e-3))

y0 <- ty[,"current"]

resODE <- try(optim(par = psi0,
                    fn = nl_try,
                    gr = gnl_try,
                    p9 = gamma.true,
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

cat(paste0("\t Fitting SDEs for MLE with bounded g_s...\n"))

res.mmfit <- try(optim(par = c(1,1), # psi0
                       fn = mm_nl,
                       gr = mm_gnl,
                       lS = log(c(7, 10.1,13.7)),
                       K0 = c(50,100,300),
                       method = "L-BFGS-B",
                       control = list(lmm = 2500,
                                      factr = 0,
                                      maxit = 10000,
                                      trace = 1,
                                      REPORT = 1),
                       hessian = FALSE),
                 silent = TRUE)

mm_sd <- sqrt(mm_s2(res.mmfit$par,
                    lS = log(c(7, 10.1,13.7)),
                    K0 = c(50,100,300)))

gs_lb <- exp(mm_fit(res.mmfit$par, 4) - qnorm(1 - 0.01/2)*mm_sd)
gs_ub <- exp(mm_fit(res.mmfit$par, 4) + qnorm(1 - 0.01/2)*mm_sd)

if(nCell == 1){
  pdf(paste0(allFigsPath, "fig11c.pdf"), width = 6, height = 5)
  par(mar = c(6,6.4,1,1.2))
  plot(cbind(c(50,100,300),
             c(7, 10.1,13.7)),
       xlim = c(0,300),
       ylim = c(0,15),
       xlab = "",
       ylab = expression(g[s]~(pS)),
       cex.axis = 2,
       cex.lab = 2,
       cex = 5,
       type="n", 
       col = 'darkgrey',
       xaxt = 'n',
       yaxt = 'n')
  lines(seq(0,500,length.out = 100),
        sapply(seq(0,500,length.out = 100),
               function(k0){
                 exp(mm_fit(res.mmfit$par, k0))
               }),
        lwd = 2)
  lines(seq(0,500,length.out = 100),
        sapply(seq(0,500,length.out = 100),
               function(k0){
                 exp(mm_fit(res.mmfit$par, k0) - qnorm(1 - 0.01/2)*mm_sd)
               }),
        lwd = 2)
  lines(seq(0,500,length.out = 100),
        sapply(seq(0,500,length.out = 100),
               function(k0){
                 exp(mm_fit(res.mmfit$par, k0) + qnorm(1 - 0.01/2)*mm_sd)
               }),
        lwd = 2)
  points(cbind(c(50,100,300),
               c(7, 10.1,13.7)),
         cex = 5,
         pch = "*",
         col = 'darkgrey')
  axis(side = 1,
       cex.axis = 2,
       cex.lab = 2,
       padj = .5)
  axis(side = 2,
       cex.axis = 2,
       cex.lab = 2,
       hadj = .5)
  title(xlab= expression("["~K~"]"[O]), line=4.5, cex.lab=2)
  dev.off()
  
  pdf(paste0(allFigsPath, "fig11c-zoom.pdf"), width = 4, height = 3.5)
  par(mar = c(3,3.1,.5,.5))
  plot(c(0.0),
       xlim = c(3.5,4.5),
       ylim = c(.8,1),
       xlab = "",
       ylab = "",
       cex.axis = 2,
       cex.lab = 2,
       cex = 5,
       type="n", 
       col = 'darkgrey',
       xaxt = 'n',
       yaxt = 'n')
  lines(seq(0,10,length.out = 100),
        sapply(seq(0,10,length.out = 100),
               function(k0){
                 exp(mm_fit(res.mmfit$par, k0))
               }),
        lwd = 2)
  lines(seq(0,10,length.out = 100),
        sapply(seq(0,10,length.out = 100),
               function(k0){
                 exp(mm_fit(res.mmfit$par, k0) - qnorm(1 - 0.01/2)*mm_sd)
               }),
        lwd = 2)
  lines(seq(0,10,length.out = 100),
        sapply(seq(0,10,length.out = 100),
               function(k0){
                 exp(mm_fit(res.mmfit$par, k0) + qnorm(1 - 0.01/2)*mm_sd)
               }),
        lwd = 2)
  points(cbind(c(4, 4),
               c(gs_lb, gs_ub)),
         cex = 5,
         pch = 20,
         lwd = 7,
         col = 'darkgrey')
  axis(side = 1,
       cex.axis = 2,
       cex.lab = 2,
       padj = .5)
  axis(side = 2,
       cex.axis = 2,
       cex.lab = 2,
       hadj = .5)
  dev.off()
}

psi0 <- c(head(resODE$par, 8),
          mm_fit(res.mmfit$par, 4) + log(1e-6),
          res.s2nu$par[1],
          resODE$par[9] - mm_fit(res.mmfit$par, 4) - log(1e-6))

res.SDE_gsfixed <- try(optim(par = psi0, # psi0
                             fn = nl.try,
                             gr = gnl.try,
                             gamma = gamma.true,
                             times = tps,
                             V = Vsin,
                             E = E.true,
                             y = y0,
                             hmax = 1, # max(tps),
                             method = "L-BFGS-B",
                             lower = c(rep(-Inf, 8),  log(gs_lb*1e-6), rep(-Inf, 2)),
                             upper = c(rep(Inf, 8),  log(gs_ub*1e-6), rep(Inf, 2)),
                             control = list(lmm = 2500,
                                            factr = 0,
                                            maxit = 10000,
                                            trace = 1,
                                            REPORT = 1),
                             hessian = FALSE),
                       silent = TRUE)

save.image(paste0(resPath,"/data-fit-nCell", nCell, ".RData"))

mus <- get.pred(psi = res.SDE_gsfixed$par,
                gamma = gamma.true,
                times = tps,
                V = Vsin,
                E = E.true,
                hmax = max(tps))

mu_psd <- (mus[1:length(Vsin) %% 10 == 1, "mu"] + qnorm(1-0.05/2)*sqrt(mus[1:length(Vsin) %% 10 == 1, "s"]))
mu_msd <- (mus[1:length(Vsin) %% 10 == 1, "mu"] - qnorm(1-0.05/2)*sqrt(mus[1:length(Vsin) %% 10 == 1, "s"]))


pdf(paste0(allFigsPath, ifelse(nCell == 5, "fig11a", paste0("fig12-cell", nCell)), ".pdf"), width = 12, height = 4.5)
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
lines(tps[1:length(tps) %% 10 == 1]/1000, 
      mus[1:length(Vsin) %% 10 == 1, "mu"],
      lwd = .5,
      col = scales::alpha("blue",.5))
legend(x = "bottomright",
       legend = c("measured current",
                  "model fit"),
       col = c("black",
               "blue"),
       lwd = 10,
       horiz = FALSE,
       bty = 'n',
       cex = 1.5)
text(x = 7.8,
     y = max(mu_psd)*.95,
     cex = 2,
     fon = 2,
     label = paste0("cell ", nCell))
dev.off()


write.csv(x = res.SDE_gsfixed$par,
          file = paste0(resPath, paste0("theta_cell_gs-fixed", nCell, ".csv")))
write.csv(x = res.SDE_gsfixed$par,
          file = paste0(resData, paste0("theta_cell_gs-fixed", nCell, ".csv")))



