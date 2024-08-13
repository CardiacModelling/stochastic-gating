
all.s2 <- c(1e-6,1e-5,1e-4)
all.nu <- c(100,1000,10000)

for (s2 in all.s2) {
  for (nu in all.nu) {
    s2_nu.true_nSim_nCores <- paste(s2, nu, 20, 20, sep = "-")

   source("./source/sim-fit.R") 
  }
}

temp = list.files(path = "./results/UQ1/", pattern="\\.csv$")
myfiles = lapply(temp, function(fn){
  read.csv(file = paste0("./results/UQ1/", fn), header = TRUE)
})


allthetas <- sapply(myfiles, function(M){t(as.matrix(M[,-1]))}, simplify = "array")
dimnames(allthetas)[[3]] <- as.vector(sapply(temp, function(s){unlist(strsplit(unlist(strsplit(s, split = "_"))[2], ".csv"))}))

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

gs.true <- g.true/nu.true

psi.true <- log(c(head(th.true, -1), gs.true, s2.true, nu.true - 1))

psitrues <- sapply(dimnames(allthetas)[[3]], function(s){
  s2_eta <- as.numeric(unlist(strsplit(s, split = "-")))
  gs <- g.true/s2_eta[2]
  psi <-  log(c(head(th.true, -1), gs, s2_eta[1], s2_eta[2] - 1))
  return(psi)
})


pdf(paste0(allFigsPath, "/fig4c.pdf"), width = 7, height = 6)
par(mar = c(5,6,2,2))
plot(c(psitrues[,rep(1:9, each = 20)]), 
     c(allthetas),
     pch = rep(c(1,2,3), each = 3*20*11),
     col = rep(rep(c("blue", "red", "green"), each = 20*11), times = 3),
     cex = 3,
     lwd = 1,
     cex.axis = 2.5,
     cex.lab = 2.5,
     xlab = expression(log~psi),
     ylab = expression(log~hat(psi)))
lines(range(c(psitrues[,rep(1:9, each = 20)],
              allthetas)),
      range(c(psitrues[,rep(1:9, each = 20)],
              allthetas)),
      lwd = 5, col = 'grey')
legend(x = -6,
       y = -11,
       legend = formatC(c(1e-6,1e-5,1e-4), format = "e", digits = 0),
       pch = 1:3,
       bty = 'n',
       horiz = TRUE,
       cex = 1.5,
       pt.cex = 3,
       pt.lwd = 3,
       text.font = 2)
legend(x = -6,
       y = -9,
       legend = formatC(c(100,1000,10000), format = "e", digits = 0),
       col = c("blue", "red", "green"),
       pch = 15,
       bty = 'n',
       horiz = TRUE,
       cex = 1.5,
       pt.cex = 4,
       pt.lwd = 3,
       text.font = 2)
dev.off()

