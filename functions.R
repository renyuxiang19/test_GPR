library(rdist)
library(tidyverse)
library(rlang)
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

# Calculate covariance matrix
# 3D
make_cov <- function(m1, m2, kernel, ... ){
  kernel <- match.fun(kernel)
  # "m1" and "m2" should be matrices or data.frames with same numbers of columns.
  #   e.g., if you want to calculate in two horizontal directions, the "m1" and "m1" should both have 2 columns.
  #   The first column in the both two matrices means one direction, and the second for another.
  ############################################################
  # Check input data.
  if (ncol(m1) != ncol(m2)){ stop("Error: The number of columns of input matrices are not consistent.") }
  # Create distance matrices
  d <- rdist::cdist(m1, m2) 
  # Calculate covariance matrices with kernel function
  k <- purrr::modify(d, kernel, ...)
  return(k)
}
