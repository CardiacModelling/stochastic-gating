
for (j in 1:11) {
  s2_nu.true_paridx_gridLength <- paste(format(0.00001, scientific = FALSE), 1000, j, 50, sep = "-")
  
  source("./source/PL-main.R")
}
