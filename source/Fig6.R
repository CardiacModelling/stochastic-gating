
for (j in c(100,
            90,
            80,
            70,
            60,
            50,
            40,
            30,
            20,
            15,
            10,
            5)) {
  s2_nu.true_nSim_nCores <- paste(format(0.00001, scientific = FALSE), 
                                  1000, 
                                  100, 
                                  20,
                                  j,
                                  sep = "-")
  
  source("./source/sim-sampFreq-main.R")
}


temp = list.files(path = respFreq, pattern="\\.csv$")
temp <- temp[order(as.numeric(sapply(sapply(sapply(temp, function(s){tail(unlist(strsplit(s, split = "-")),1)}), 
       function(s){unlist(strsplit(s, split = "pFreq"))[2]}), function(s){unlist(strsplit(s, split = ".csv"))})))]
myfiles = lapply(temp, function(fn){
  read.csv(file = paste0(respFreq, fn), header = TRUE)
})


all.estimates <- sapply(myfiles, function(myfile){as.matrix(myfile[,-1])}, simplify = "array")
dimnames(all.estimates)[[3]] <- c(5,10,15,20,30,40,50,60,70,80,90,100)

pal.col <- c("#A6EAD6",
             "#EA8CB6",
             "#E4EA8B",
             "#A5CEE3",
             "#E2E6DF",
             "#EEB691",
             "#DAE5B6",
             "#D789EB",
             "#A8EC9B",
             "#CFB0AD",
             "#AFA7E4",
             "#E7C2E5")

pdf(paste0(allFigsPath, "/fig6a.pdf"), width = 6, height = 5)
par(mar = c(4.5,5.5,.5,.5))
plot(rep(psi.true,100), 
     c(t(all.estimates[,,1])),
     cex = 2,
     lwd = 3,
     col = scales::alpha(pal.col[1],.2),
     xlab = expression(psi[true]),
     ylab = expression(hat(psi)),
     cex.axis = 2,
     cex.lab = 2,
     xlim = range(all.estimates, psi.true),
     ylim = range(all.estimates, psi.true))
sapply(2:12, function(pfreq){
  points(rep(psi.true,100), 
         c(t(all.estimates[,,pfreq])),
         cex = 2,
         lwd = 3,
         col = scales::alpha(pal.col[pfreq],.2))
})
lines(range(all.estimates, psi.true),
      range(all.estimates, psi.true),
      col = scales::alpha("grey",.5),
      lwd = 5)
legend(x = "bottomright",
       legend = paste0(dimnames(all.estimates)[[3]],"%"),
       col = pal.col,
       bty = 'n',
       pch = 1,
       cex = 1.2,
       pt.cex = 2,
       pt.lwd = 3)
dev.off()


pdf(paste0(allFigsPath, "/fig6b.pdf"), width = 15, height = 5)
par(mar = c(2,5.5,.5,.5))
boxplot(sapply(1:dim(all.estimates)[3], function(k){
  apply(all.estimates[,,k], 1, FUN = function(x){
    norm(x - psi.true, type = "2")/norm(psi.true, type = "2")
  })
}), 
outline = FALSE,
col = scales::alpha(pal.col,.7),
lwd = 3,
cex.axis = 2,
cex.lab = 2,
xaxt = 'n')
axis(side = 1,
     at = 1:12,
     labels = paste0(dimnames(all.estimates)[[3]],"%"),
     cex.axis = 2,
     cex.lab = 2)
dev.off()


