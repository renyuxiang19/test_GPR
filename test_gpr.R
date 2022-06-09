library(tidyverse)
library(rdist)
library(rSPDE)
# Training data
s <- c(0.1, 0.3, 0.5, 0.7, 0.9)
z <- c(6, 4, 4, 7, 6)
z<- as.matrix(z)
par(mar=c(5,4,1,1)+0.1)
plot(s,z, 
     ylim = c(0,8), pch=19)

# Kernel functions
sof <- 0.5
sd <- 5
nu <- 0.5
noise <- TRUE

kernel_g <- function(d){
  k <- (sd^2) * exp(-pi* (d/sof)^2) 
  return(k)
}

kernel_m <- function(d){
  k <- (sd^2) * exp(-2* (d/sof))
  return(k)
}

kernel_b <- function(d){
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

kernel_wm <- function(d){
  if (d == 0) {
    covariance <- sd^2
  }else{
    frac <- sqrt(pi)*gamma(nu+0.5)*d/(gamma(nu)*sof)
    covariance <-sd^2 * (2/gamma(nu)) * frac^nu * besselK(2*frac, nu)
  }
  return(covariance)
}
kernel_wm2 <- function(d,k){
  if (d == 0) {
    covariance <- sd
  }else{
    covariance <-sd^2 * (2^(1-nu)/gamma(nu)) * (k*d)^nu * besselK(k*d, nu)
  }
  return(covariance)
}
# Calculate covariance matrices
make_cov <- function(s1, s2, kernel){
  distances <- rdist::cdist(s1,s2,metric = "manhattan")
  cov_m <- switch (kernel,
    Gaussian = kernel_g(distances),
    Markovian = kernel_m(distances),
    Binary = kernel_b(distances),
    Whittle = purrr::modify(distances, 
                            rSPDE::matern.covariance, 
                            kappa = 10, nu = nu, sigma = sd),
    wm2 = purrr::modify(distances, kernel_wm)
  )
  return(cov_m)
}

# Predict
test_s <- seq(0, 1, length= 40)

k <- make_cov(s1=s, s2=s, kernel = "Markovian")
k_star <- make_cov(s1=test_s, s2=s, kernel = "Markovian")
 # Add noise 
if (noise){
  r <- diag(nrow(k))
  k <- k + r
}


z_predict <- k_star %*% solve(k) %*% z

# Plot predicting
lines(test_s,z_predict)

