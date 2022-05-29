# Test GPR with simple 3D data.
library(tidyverse)
library(rdist)

# training data
v1 <- c(1,3,4,6,7)
v2 <- c(2,4,6,5,3)
x1 <- rep(1,length(v1))
x2 <- rep(5,length(v1))
z1 <- seq_along(v1)
z2 <- seq_along(v2)
training <- data.frame(x=c(x1,x2), z=c(z1,z2), value=c(v1,v2))
training$y <- 0
training[c("x","y")]
testing <- data.frame(x = rep(1:5, each=40), z = rep(seq(1,5,length=40), times=5)
  , y = rep(0,5*40))
#
train_pic <- ggplot(data=training) +
        geom_point(mapping = aes(x=value,y=z))+
        geom_line(mapping = aes(x=value,y=z),color="#D3323F",orientation = "y")+
        scale_y_reverse()+
        facet_wrap(~x)
train_pic
#
make_cov_3D <- function(m1, m2, sof_h=5, sof_v=5, sd0 = 5){
  # Create distance matrices
  d_h <- rdist::cdist(m1[c("x","y")],m2[c("x","y")]) 
  d_v <- rdist::cdist(m1$z,m2$z)
  # Calculate covariance matrices in two directions.
  k_h <- modify(d_h, kernel_g, sof = sof_h, sd = sd0)
  k_v <- modify(d_v, kernel_g, sof = sof_v, sd = sd0)
  # Calculate the 3D covariance matrix
  k <- k_h * k_v
  return(k)
}
#

k11 <- make_cov_3D(m1 = training, m2 = training)
k21 <- make_cov_3D(m1 = testing, m2 = training)

noise <- rnorm(length(k21),mean=0,sd=1) %>% matrix(ncol = ncol(k21), nrow = nrow(k21)) 
r <- cov(noise)

testing$value <-  k21 %*% solve(k11+r) %*% training$value

test_pic <- ggplot() +
        geom_point(data=training,mapping = aes(x=value,y=z))+
        geom_line(data=testing, mapping = aes(x=value,y=z),color="#D3323F",orientation = "y")+
        scale_y_reverse()+
        facet_wrap(~x)
test_pic
