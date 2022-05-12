library(tidyverse)
source("functions.r")
#
n_sws <- read_MAIC("kaminokoike_SWS.dat")

x <- as.numeric(attr(factor(n_sws$X),"levels")) 
training <- list()
training$z <- split(n_sws$Z, n_sws$X)
training$nsws <- split(n_sws$Nsws, n_sws$X)
predi <- list()
for (i in seq_along(x)) {
  predi[[i]] <- gpr_1d(x = training$z[[i]], y = training$nsws[[i]], kernel = kernel_g, sof = 3, sd0 = 5)
}

i <- 8
par(mar=c(5,4,1,1)+0.1)
plot(training$nsws[[i]],training$z[[i]], 
     ylim = c(0,8), pch=19, type = "l", lty = 1)
lines(predi[[i]]$y, predi[[i]]$x, lty=2)
