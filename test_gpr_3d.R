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
# 
s1 <- c(z1, 10 + z2)
s2 <- 1:15

make_cov <- function(shi1, shi2, sv1, sv2, kernel, ...){
  len_s1 <- length(sh1)
  len_s2 <- length(sh2)
  cov_m <- matrix(nrow = len_s2, ncol = len_s1)
  for (i in 1:len_s1) {
    cov_m[,i] <- map_dbl(s1, function(x) kernel(s1=x, s2=s2[i], ...) )
  }
  return(cov_m)
}
make_cov <- functon
#
v <- matrix(c(v1,v2),ncol = 2)
k <- make_cov(s1=x, s2=x, kernel = kernel_wm, ...)
b <- data.frame(z1,z2)
a <- rdist::cdist(v,v)
dim(a)
sloop::s3_class(v1)
class(a)
attributes(a)

switch(3, 10,20,30,40,50)
switch("cc", a = 1, cc =, cd =3, d = 2)
