library(rdist)
library(tidyverse)
# Read_data
read_MAIC <- function(file_name){
  rawdat<-read.csv(file = file_name) #カンマ区切りのファイルを読み込む(他の書式はMAICの入力fileと同じ)　
  colnames(rawdat)<-c("Z","Nsws")
  pointnum<-rawdat[1,1]
  indexnum<-numeric(length = pointnum)
  xcoordinate<-numeric(length = pointnum)
  indexnum[1]<-rawdat[2,1]
  xcoordinate[1]<-rawdat[3,1]
  maindata<-rawdat[4:(3+indexnum[1]),]
  cumdatarow<-0
  for (i in 2:pointnum) {
    cumindexrow<-1+(i-1)*2
    cumdatarow<-indexnum[(i-1)]+cumdatarow
    indexnum[i]<-rawdat[cumdatarow+cumindexrow+1,1]
    xcoordinate[i]<-rawdat[cumdatarow+cumindexrow+2,1]
    maindata<-rbind(maindata,rawdat[(cumdatarow+cumindexrow+3):(cumdatarow+cumindexrow+2+indexnum[i]),]) 
  }
  if(length(rep(xcoordinate,indexnum))!=nrow(maindata)){
    print("Error")
  }else{maindata$X<-rep(xcoordinate,indexnum)}
  return(maindata)
}

# Kernel functions
kernel_g <- function(s1, s2, sof, sd0){
  k <- (sd0^2) * exp(-pi* (abs(s1-s2)/sof)^2) 
  return(k)
}

kernel_m <- function(s1, s2,sof, sd){
  k <- (sd0^2) * exp(-2* (abs(s1-s2)/sof))
  return(k)
}

kernel_b <- function(s1, s2, sof, sd0){
  dif <- abs(s1-s2)
  if (dif < sof) {
    k <- (sd0^2) * (1-(dif/sof)) 
  }else{
    k <- 0
  }
  return(k)
}

kernel_e <- function(s1, s2){
  k <- exp(-0.5*abs(s1-s2)^2)
  return(k)
}
# 3D Whittle-Matern kernel 
d_h <- function(shi1, shi2, shj1, shj2){
  sqrt((shi1-shi2)^2+(shj1-shj2)^2)
}

d_v <- function(sv1,sv2){
  abs(sv1-sv2)
}

k_w <- function(d, nu, sof){
  frac <- sqrt(pi)*gamma(nu+0.5)*d/(gamma(nu)*sof)
  ans <- (2/gamma(nu)) * frac^nu * besselK(2*frac, nu)
  return(ans)
}

kernel_wm <- function(shi1,shi2,shj1=0,shj2=0,sv1,sv2, nu_h, sof_h, nu_v, sof_v){
  dh <- d_h(shi1, shi2, shj1, shj2)
  dv <- d_v(sv1, sv2)
  ans <- k_w(dh, nu_h, sof_h) * k_w(dv, nu_v, sof_v)
  return(ans)
}

# Calculate covariance matrices
make_cov <- function(s1, s2, kernel, ...){
  len_s1 <- length(s1)
  len_s2 <- length(s2)
  cov_m <- matrix(nrow = len_s1, ncol = len_s2)
  for (i in 1:len_s2) {
    cov_m[,i] <- map_dbl(s1, function(x) kernel(s1=x, s2=s2[i], ...) )
  }
  return(cov_m)
}

# GPR functions 
gpr_1d <- function(x, y, num=100, kernel, sd_noise=1, ...){
  # num : number of test data
  # x, y: training data, which should be vectors with same length.
  # kernel : kernel function.
  if(length(x)!=length(y)){
    stop("The lengths of 'x' and 'y' should be same")
  }else{
    len_tra <- length(x)
  }
  test_x <- seq(min(x), max(x), length = num)
  # Calculate covariance matrices
  k <- make_cov(s1=x, s2=x, kernel = kernel, ...)
  k_star <- make_cov(s1=test_x, s2=x, kernel = kernel, ...)
  # Add noise 
  noise <- rnorm(length(k_star),mean=0, sd=sd_noise) %>% matrix(ncol = ncol(k_star), nrow = nrow(k_star))
  r <- cov(noise)
  # predict y
  predict_y <- k_star %*% solve(k+r) %*% y
  ans <- data.frame(x = test_x, y = predict_y)
  return(ans)
}

