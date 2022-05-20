library(tidyverse)
library(rdist)
library(sloop)
source("functions.r")
# training data
v1 <- c(1,3,4,6,7)
v2 <- c(2,4,6,5,4)
x1 <- rep(1,length(v1))
x2 <- rep(3,length(v1))
z1 <- seq_along(v1)
z2 <- seq_along(v2)
training <- data.frame(x=c(x1,x2), z=c(z1,z2), value=c(v1,v2))
training$y <- 0
training[c("x","y")]
# 

make_cov <- function(training, nu_h=0.5, nu_v=0.5, sof_h=1, sof_v=1, sd = 5){
 d_h <- rdist::cdist(training[c("x","y")],training[c("x","y")]) 
 d_v <- rdist::cdist(training$z,training$z)
 k_h <- kernel_wm(d_h, nu_h, sof_h, sd)
 k_v <- kernel_wm(d_v, nu_v, sof_v, sd)
 browser()
}

make_cov(training)
#
