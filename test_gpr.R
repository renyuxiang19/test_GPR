library(tidyverse)
library(rdist)
# Training data
s <- c(0.1, 0.3, 0.5, 0.7, 0.9)
z <- c(6, 4, 4, 7, 6)
z<- as.matrix(z)
par(mar=c(5,4,1,1)+0.1)
plot(s,z, 
     ylim = c(0,8), pch=19)

# Kernel functions
sof <- 5
sd <- 5

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

# Calculate covariance matrices
make_cov <- function(s1, s2, kernel){
  distances <- rdist::cdist(s2,s1,metric = "manhattan")
  cov_m <- switch (kernel,
    Gaussian = kernel_g(distances),
    Markovian = kernel_m(distances),
    Binary = kernel_b(distances)
  )
  return(cov_m)
}

# Predict
test_s <- seq(0, 1, length= 40)

k <- make_cov(s1=s, s2=s, kernel = "Gaussian")
k_star <- make_cov(s1=s, s2=test_s, kernel = "Gaussian")
 # Add noise 
noise <- rnorm(length(k_star),mean=0,sd=1) %>% matrix(ncol = ncol(k_star), nrow = nrow(k_star))
r <- cov(noise)

z_predict <- k_star %*% solve(k+r) %*% z

# Plot predicting
lines(test_s,z_predict)
