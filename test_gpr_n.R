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
  predi[[i]] <- gpr_1d(x = training$z[[i]], y = training$nsws[[i]], kernel = kernel_g, sof = 10, sd0 = 20)
}

i <- 4
par(mar=c(5,4,1,1)+0.1)
plot(training$nsws[[i]],training$z[[i]],
     ylim = rev(range(training$z[[i]])), pch=19, type = "l", lty = 1,)
lines(predi[[i]]$y, predi[[i]]$x, lty=2)

# meta_data <- transpose(training)
# 
# ggplot() +
#   geom_point(data = data.frame(meta_data[[i]]), mapping = aes(x=nsws,y=z))+
#   geom_line(data = data.frame(predi[[i]]), mapping = aes(x=y,y=x),color="#D3323F")+
#   scale_y_reverse()+
#   xlim(0,4)
# 
# data.frame(predi[[i]])
