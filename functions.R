library(rdist)
library(tidyverse)
library(rlang)
# Read_data
read_MAIC <- function(file_name){
  rawdat<-read.csv(file = file_name) #カンマ区切りのファイルを読み込む(他の書式はMAICの入力fileと同じ)　
  colnames(rawdat)<-c("z","nsws")
  point_num<-rawdat[1,1]
  #
  data_number_loca <- 2
  data_number <- numeric(length = point_num)
  distance <- numeric(length = point_num)
  for (i in 1:point_num){
    data_number[i] <- rawdat[data_number_loca, 1]
    distance[i] <- rawdat[data_number_loca + 1 , 1]
    data_number_loca <- data_number_loca + data_number[i] + 2
    # check the number of data.
    if (!is_mathinteger(data_number[i])) {
      paste("Error: the number of data at ", i,
            "-th point is not a integer, pleas check the number of data in input file.",sep = "") |>
        stop()
    }
  }
  distance_index <- purrr::map2(distance, data_number, function(x, y) c(NA, NA, rep(x, each = y)) ) |> 
    unlist()
  distance_index <- c(NA, distance_index)
  #
  rawdat <- rawdat |> dplyr::mutate(x = distance_index)
  dplyr::filter(rawdat, !is.na(x))
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

# Whittle-Matern kernel 
kernel_wm <- function(d, nu, sof, sd){
  if (d == 0) {
    covariance <- sd^2
  }else{
    frac <- sqrt(pi)*gamma(nu+0.5)*d/(gamma(nu)*sof)
    covariance <- sd^2 * (2/gamma(nu)) * frac^nu * besselK(2*frac, nu)
  }
  return(covariance)
}

# Calculate covariance matrix
make_cov <- function(m1, m2, kernel, ... ){
  # "m1" and "m2" should be matrices or data.frames with same numbers of columns.
  #   e.g., if you want to calculate in two horizontal directions, the "m1" and "m1" should both have 2 columns.
  #   The first column in the both two matrices means one direction, and the second for another.
  ############################################################
  kernel <- match.fun(kernel)
  # Check input data.
  if (ncol(m1) != ncol(m2)){ stop("Error: The number of columns of input matrices are not consistent.") }
  # Create distance matrices
  d <- rdist::cdist(m1, m2) 
  # Calculate covariance matrices with kernel function
  k <- purrr::modify(d, kernel, ...)
  return(k)
}
#
is_mathinteger = function(v){
  all(round(v) == v)
}

