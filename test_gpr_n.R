# Test GPR with real N-value data.
# library(tidyverse)
# library(rlang)
library(MASS)
source("functions.r")

## Parameters of GPR
filename <- "kaminokoike_SWS.dat"
kernel_fun <- "kernel_g"
sof_h <- 10
sof_v <- 6
sd <- 2
nu <- 1.25
noise <- TRUE

### Unimportant parameter
mesh_size_v <- 0.1

## read 2D data.
n_sws <- read_MAIC(filename) |> mutate(y=0)

## Plot raw data.
raw_pic <- ggplot(data=n_sws)+
  geom_point(mapping = aes(x=nsws,y=z))+
  geom_line(mapping = aes(x=nsws,y=z),orientation = "y")+
  scale_y_reverse()+
  facet_wrap(~x, nrow = 2)

## Prepare testing data (mesh).
depth <- n_sws |> 
  group_by(x) |> 
  summarise(min_depth = min(z),max_depth = max(z))
testing <- list(depth$x, depth$min_depth, depth$max_depth) |> 
  pmap_dfr(function(dis, min, max){
    value <- seq(min, max, by = mesh_size_v)
    len <- length(value)
    tibble(x = rep(dis,times = len), z = value)
    }) |>
  mutate(y=0)

## Calculate cov matrix..
make_k11 <- function(noise = TRUE){
  ## Calculate cov matrix.
  k11_h <- make_cov(m1 = n_sws[c("x","y")], m2 = n_sws[c("x","y")],
                    kernel = kernel_fun, sof = sof_h, sd = sd)
  k11_v <- make_cov(m1 = n_sws["z"], m2 = n_sws["z"],
                    kernel = kernel_fun, sof = sof_v, sd = sd)
  k11 <- k11_h * k11_v
  if (noise){
    r <- diag(nrow(k11))
    k11 <- k11 + r
  }
  return(k11)
}
make_k21 <- function(){
  k21_h <- make_cov(m1 = testing[c("x","y")], m2 = n_sws[c("x","y")],
                    kernel = kernel_fun, sof = sof_h, sd = sd)
  k21_v <- make_cov(m1 = testing["z"], m2 = n_sws["z"],
                    kernel = kernel_fun, sof = sof_v, sd = sd)
  k21 <- k21_h * k21_v
  return(k21)
}
k21 <- make_k21()
k11 <- make_k11( noise = noise )

## Predict
testing$nsws <- k21 %*% ginv(k11) %*% n_sws$nsws

## plot predicting
test_pic <- ggplot() +
  geom_point(data=n_sws,mapping = aes(x=nsws,y=z))+
  geom_line(data=testing, mapping = aes(x=nsws,y=z),color="#D3323F",orientation = "y")+
  scale_y_reverse()+
  facet_wrap(~x, nrow = 2)
test_pic
