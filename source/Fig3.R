cat("\nInstall/load packages")
inst.pkgs <- installed.packages()
l.pkgs <- c("Deriv",
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
source("./source/sim-example-Fig3-fun.R")

set.seed(123)

s2_nu <- "0-100"
# s2_nu <- commandArgs(trailingOnly = TRUE) # s2_nu <- "0-100"
s2.true <- as.numeric(unlist(strsplit(s2_nu, split = "-"))[1])
nu.true <- as.numeric(unlist(strsplit(s2_nu, split = "-"))[2])

##############
## currTime ##
##############

currTime <- Sys.time()
currTime <- gsub(" ", "_", currTime)
currTime <- gsub("-", "", currTime)
currTime <- gsub(":", "", currTime)
currTime <- paste(currTime, s2_nu, sep = "_")

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

# gamma0 <- 1/6.7
# delta0 <- 1/2.5

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

gamma.true <- tail(th.true, 1)


Y0 <- get.x.em(theta = th.true,
               nu = nu.true,
               times = tps,
               chs.lst = chs)

X0 <- get.logistic(Y0)

y0 <- get.y.em(times = tps,
               X = X0,
               gs = g.true/nu.true,
               nu = nu.true,
               E = E.true,
               s2 = s2.true)

psi.true <- log(c(head(th.true, -1), 
                  g.true/nu.true, 
                  0, 
                  nu.true - 1))

mus <- get.pred(psi = psi.true,
                gamma = gamma.true,
                times = tps,
                V = Vsin,
                E = E.true,
                hmax = max(tps))

X.pred <- get.pred.X(psi = c(head(psi.true,-1), 0),
                     gamma = gamma.true,
                     times = tps,
                     V = Vsin,
                     E = E.true,
                     hmax = 1)

VAR_X <- X.pred[1:nrow(X.pred) %% 10 == 1, c(1,5,8,10)+4]
COV_X <- X.pred[1:nrow(X.pred) %% 10 == 1, c(2,3,4,6,7,9)+4]
VAR_X_last <- rowSums(cbind(2*COV_X, VAR_X))

X_msd <- (cbind(X.pred[1:nrow(X.pred) %% 10 == 1,1:4], 1 - rowSums(X.pred[1:nrow(X.pred) %% 10 == 1,1:4])) 
        - qnorm(1 - 0.05/2)*sqrt(cbind(VAR_X, VAR_X_last)*diff(tps)[1]/nu.true))
X_psd <- (cbind(X.pred[1:nrow(X.pred) %% 10 == 1,1:4], 1 - rowSums(X.pred[1:nrow(X.pred) %% 10 == 1,1:4])) 
        + qnorm(1 - 0.05/2)*sqrt(cbind(VAR_X, VAR_X_last)*diff(tps)[1]/nu.true))

mu_psd <- (mus[1:length(Vsin) %% 10 == 1, "mu"] + qnorm(1 - 0.05/2)*sqrt(mus[1:length(Vsin) %% 10 == 1, "s"]))
mu_msd <- (mus[1:length(Vsin) %% 10 == 1, "mu"] - qnorm(1 - 0.05/2)*sqrt(mus[1:length(Vsin) %% 10 == 1, "s"]))

pdf(paste0(allFigsPath, "fig3.pdf"), width = 12, height = 8)
layout(mat = matrix(data = c(1,2,2,3,3), ncol = 1))
par(mar = c(1,5,1,1))
plot(tps[1:length(tps) %% 10 == 1]/1000,
     Vsin[1:length(Vsin) %% 10 == 1],
     type = 'l', 
     lty = 1,
     lwd = 1,
     xlab = "",
     ylab = "V (mV)",
     cex.axis = 2,
     cex.lab = 2,
     xaxt = 'n')
###
plot(1, 
     type="n", 
     xlab="", 
     ylab="x", 
     xlim=c(0, 8), 
     ylim=c(0, 1),
     cex.axis = 2,
     cex.lab = 2,
     col = scales::alpha(1:5,.5),
     xaxt = 'n')
sapply(1:5, function(j){
  polygon(c(tps[1:length(tps) %% 10 == 1]/1000, 
          rev(tps[1:length(tps) %% 10 == 1]/1000)), 
        c(X_psd[,j], rev(X_msd[,j])),
        col = scales::alpha(j,.2), lty = 0)
})
###
matplot(tps[1:length(tps) %% 10 == 1]/1000,
        cbind(X0[1:nrow(X0) %% 10 == 1,], 1 - rowSums(X0[1:nrow(X0) %% 10 == 1,])), 
        type = 'l', 
        lty = 1,
        lwd = .5,
        xlab = "",
        ylab = "x",
        cex.axis = 2,
        cex.lab = 2,
        col = scales::alpha(1:5,.5),
        xaxt = 'n',
        add = TRUE)
###
matplot(tps[1:length(tps) %% 10 == 1]/1000,
        cbind(X.pred[1:nrow(X.pred) %% 10 == 1,1:4], 1 - rowSums(X.pred[1:nrow(X.pred) %% 10 == 1,1:4])), 
        type = 'l', 
        lty = 1,
        lwd = 1,
        cex.axis = 2,
        cex.lab = 2, 
        add = TRUE)
legend(x = "topright",
       legend = chs,
       col = 1:5,
       lwd = 10,
       horiz = TRUE,
       bty = 'n',
       cex = 1.5)
legend(x = 7,
       y = .3,
       legend = c("simulation",
                  "DMEs"),
       col = c("black",
               "black"),
       lwd = c(1,5),
       horiz = FALSE,
       bty = 'n',
       cex = 1.8)
###
par(mar = c(2,5,1,1))
###
plot(1, 
     type="n", 
     xlab = "time (s)",
     ylab = "y (nA)",
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
     y0[1:length(y0) %% 10 == 1], 
     lty = 1,
     lwd = .5,
     xlab = "time (s)",
     ylab = "y (nA)",
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


