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
source("./source/PL-fun.R")

# s2_nu.true_paridx_gridLength <- commandArgs(trailingOnly = TRUE) # s2_nu.true_paridx_gridLength <- "0.00001-1000-1-50"
s2.true <- as.numeric(unlist(strsplit(s2_nu.true_paridx_gridLength, split = "-"))[1])
nu.true <- as.numeric(unlist(strsplit(s2_nu.true_paridx_gridLength, split = "-"))[2])
par.idx <- as.numeric(unlist(strsplit(s2_nu.true_paridx_gridLength, split = "-"))[3])
gridLength <- as.numeric(unlist(strsplit(s2_nu.true_paridx_gridLength, split = "-"))[4])

cat(paste0("\ts2 = ", s2.true, 
           "\tnu = ", nu.true, 
           "\tnpar.idx = ", par.idx, 
           "\tgridLength = ", gridLength, "\n"))

##############
## currTime ##
##############

# currTime <- print(as.numeric(Sys.time())*1000, digits=20)
currTime <- Sys.time()
currTime <- gsub(" ", "_", currTime)
currTime <- gsub("-", "", currTime)
currTime <- gsub(":", "", currTime)
currTime <- paste(currTime, s2_nu.true_paridx_gridLength, sep = "_")

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
                   # theta = th.true,
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
                     hessian = TRUE),
               silent = TRUE)
H_mle <- res.SDE$hessian

psi.ranges <- cbind(res.SDE$par - res.SDE$par/10, 
                    res.SDE$par + res.SDE$par/10)
psi.ranges <- t(apply(psi.ranges, 1, sort, decreasing = FALSE))

psi.range.idx <- psi.ranges[par.idx,]
I_left <- rev(head(seq(from = psi.range.idx[1],
                       to = res.SDE$par[par.idx],
                       length.out = gridLength + 1), -1))

I_right <- seq(from = res.SDE$par[par.idx], 
               to = psi.range.idx[2], 
               length.out = gridLength + 1)[-1]


PL_left <- get.PL(phi = res.SDE$par,
                  i = par.idx,
                  IC_i = I_left,
                  gamma = gamma.true,
                  times = tps,
                  V = Vsin,
                  E = E.true,
                  y = y0,
                  hmax = 1)

PL_right <- get.PL(phi = res.SDE$par,
                   i = par.idx,
                   IC_i = I_right,
                   gamma = gamma.true,
                   times = tps,
                   V = Vsin,
                   E = E.true,
                   y = y0,
                   hmax = 1)

save.image(paste0(resPath,"/PL-fit-", par.idx, ".RData"))


PL_all <- rbind(PL_left$PL.values[length(I_left):1, c(par.idx, ncol(PL_left$PL.values))],
                c(res.SDE$par[par.idx], res.SDE$value),
                PL_right$PL.values[, c(par.idx, ncol(PL_left$PL.values))])

write.csv(x = PL_all,
          file = paste0(resPath, "PL_values-", par.idx, ".csv"))


fim_mle <- get.fim(psi = res.SDE$par,
                   gamma = gamma.true,
                   times = tps,
                   V = Vsin,
                   E = E.true,
                   hmax = 1)

PL_all_complete <- rbind(PL_left$PL.values[length(I_left):1, 1:length(psi.true)],
                         c(res.SDE$par),
                         PL_right$PL.values[, 1:length(psi.true)])

fim_x <- sort(c(seq(min(I_left),
                    max(I_right),
                    length.out = 100), psi.true[par.idx]))
fim_par <- sapply(1:length(fim_x), function(j){
  psi_i <- rep(NA, length(res.SDE$par))
  
  psi_i[-par.idx] <- res.SDE$par[-par.idx]
  psi_i[par.idx] <- fim_x[j]
  
  return(FIM_parabola(x = psi_i, 
                      a = res.SDE$par, 
                      H = fim_mle, # H_mle,
                      gamma = gamma.true,
                      times = tps,
                      V = Vsin,
                      E = E.true,
                      y = y0,
                      hmax = 1))
})

save.image(paste0(resPath,"/PL-fit-", par.idx,".RData"))

write.csv(x = cbind(fim_x, fim_par),
          file = paste0(resPath, "FIM_values-", par.idx,".csv"))

pdf(paste0(allFigsPath, "PL_psi", par.idx,".pdf"), width = 7, height = 4.5)
par(mar = c(3.5,4,1,1))
plot(PL_all[,1],
     PL_all[,2]/10000,
     type = 'l',
     lwd = 10,
     ylim = c(-76.17258, -71.54193),
     # pch = 20,
     cex.axis = 3,
     cex.lab = 3,
     xlab = "",
     ylab = "",
     xaxt = 'n',
     yaxt = 'n',
     bty = 'n')
axis(side = 1,
     cex.axis = 3,
     padj = .5)
axis(side = 2,
     cex.axis = 3,
     hadj = 0)
lines(fim_x,
      fim_par/10000,
      lwd = 10,
      col = "grey")
par(xpd=TRUE)
points(res.SDE$par[par.idx],
       res.SDE$value/10000,
       pch = "*",
       col = "blue",
       cex = 10)
dev.off()
