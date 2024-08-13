
all.s2 <- c(1e-6,1e-5,1e-4)
for (s2 in all.s2) {
  s2_nu.true.min.max_nSim_nCores <- paste(format(s2, scientific = FALSE), 100, 10000, 100, 10, sep = "-")
  source("./source/sim-fit-nCh.R")
}

  