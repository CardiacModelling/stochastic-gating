
for (nCell in 1:9) {
  source("./source/data-fit-main.R")
}

thetas <- sapply(1:9, function(nCell){
  theta <- read.csv(file = paste0(resData, paste0("theta_cell_gs-fixed", nCell, ".csv")))[,2]

  return(theta)
}, simplify = "array")


x <- lapply(1:11, function(j){thetas[j,]})

pdf(paste0(allFigsPath, "fig11b.pdf"), width = 12, height = 4)
par(mar = c(3.2,4,.5,.5))
plot(0, 
     0, 
     type="n", 
     ylim=range(thetas), 
     xlim=c(1, nrow(thetas)), 
     xaxt = 'n', 
     xlab ="", 
     ylab = "",
     main ="",
     cex.axis = 3,
     cex.lab = 2)

sapply(1:nrow(thetas), function(j){
  sapply(1:ncol(thetas), 
         function(nc){
           stripchart(x[[j]][nc],
                      method="jitter",
                      pch=16,
                      col = scales::alpha(rainbow(ncol(thetas))[nc], .5),
                      add = T,
                      at = j,
                      cex = 2,
                      vertical = TRUE)})
})
axis(side = 1,
     at = 1:nrow(thetas),
     labels = head(c(sapply(1:8, function(j){bquote(theta[.(j)])}), 
                     expression(g[s]),
                     expression(sigma^2),
                     expression(eta),
                     expression(1)), -1),
     padj = .5,
     cex.axis = 3)
legend(x = "topleft",
       horiz = TRUE,
       legend = 1:9,
       col = scales::alpha(rainbow(9), .5),
       lwd = 10,
       cex = 1.6,
       bty = 'n')
dev.off()

xtable::xtable(x = rbind(thetas,
                         exp(thetas[9,])*(1 + exp(thetas[11,]))))

write.csv(x = thetas,
          file = paste0(resData, paste0("table3.csv")))
