library(rdist)
library(tidyverse)
# Read_data
read_MAIC <- function(file_name){
  rawdat<-read.csv(file = file_name) #カンマ区切りのファイルを読み込む(他の書式はMAICの入力fileと同じ)　
  colnames(rawdat)<-c("z","nsws")
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
  }else{maindata$x<-rep(xcoordinate,indexnum)}
  return(maindata)
}

# Kernel functions
kernel_g <- function(d,sof,sd){
  k <- (sd^2) * exp(-pi* (d/sof)^2) 
  return(k)
}

kernel_m <- function(d,sof,sd){
  k <- (sd^2) * exp(-2* (d/sof))
  return(k)
}

kernel_b <- function(d,sof,sd){
  dif <- d
  if (dif < sof) {
    k <- (sd^2) * (1-(dif/sof)) 
  }else{
    k <- 0
  }
  return(k)
}

kernel_e <- function(d){
  k <- exp(-0.5*d^2)
  return(k)
}

# 3D Whittle-Matern kernel 
kernel_wm <- function(d, nu, sof, sd){
  if (d == 0) {
    covariance <- sd
  }else{
    frac <- sqrt(pi)*gamma(nu+0.5)*d/(gamma(nu)*sof)
    covariance <- (2/gamma(nu)) * frac^nu * besselK(2*frac, nu)
  }
  return(covariance)
}
# Calculate covariance matrices
# make_cov <- function(s1, s2, kernel, ...){
#   len_s1 <- length(s1)
#   len_s2 <- length(s2)
#   cov_m <- matrix(nrow = len_s1, ncol = len_s2)
#   for (i in 1:len_s2) {
#     cov_m[,i] <- map_dbl(s1, function(x) kernel(s1=x, s2=s2[i], ...) )
#   }
#   return(cov_m)
# }
# 3D
make_cov <- function(m1, m2, kernel, ... ){
  # "m1" and "m2" should be matrices or data.frames with same number of columns.
  #   If you want to calculate two horizontal directions, the "m1" and "m1" should both have 2 columns.
  #   The first column in the two matrices means one direction, and the second for another.
  ############################################################
  # Check input data.
  if (ncol(m1) != ncol(m2)){ stop("Error: The number of columns of input matrices are not consistent.") }
  # Create distance matrices
  d <- rdist::cdist(m1, m2) 
  # Calculate covariance matrices with kernel function
  k <- purrr::modify(d, kernel, ...)
  return(k)
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

