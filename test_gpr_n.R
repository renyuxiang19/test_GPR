# Test GPR with real N-value data.
# library(tidyverse)
# library(rlang)
library(MASS)
library(optimx)
source("functions.r")

## Parameters of GPR
filename <- "kaminokoike_SWS.dat"
kernel_fun <- "kernel_m"
sof_h_t <- 5
sof_v_t <- 8
sd_t <- 5
sof_h_r <- 10
sof_v_r <- 10
sd_r <- 10
nu <- 1.25
noise <- TRUE

### Unimportant parameter
mesh_size_v <- 0.1

## read 2D data.
n_sws <- read_MAIC(filename) |> 
  mutate(y=0) |>
  rowid_to_column("ID") |>
  select(-ID)

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

## Calculate cov matrix.
### Using log-likelihood to optimize parameters
#### Observation vector.
z <- n_sws$nsws |> matrix()
#### Parameter vector.
para <- c(sof_h_t, sof_v_t, sd_t, sof_h_r, sof_v_r, sd_r)
#### K11 function.
make_k11 <- function(para){
  sof_h <- para[1]
  sof_v <- para[2]
  sd <- para[3]
  k11_h <- make_cov(m1 = n_sws[c("x","y")], m2 = n_sws[c("x","y")],
                    kernel = kernel_fun, sof = sof_h, sd = sd)
  k11_v <- make_cov(m1 = n_sws["z"], m2 = n_sws["z"],
                    kernel = kernel_fun, sof = sof_v, sd = sd)
  k11 <- k11_h * k11_v
  return(k11)
}
#### likelihood function.
ln_likelihood <- function(para){
  sof_h_t <- para[1]
  sof_v_t <- para[2]
  sd_t <- para[3]
  sof_h_r <- para[4]
  sof_v_r <- para[5]
  sd_r <- para[6]
  #
  m <- length(para)
  k11_t <- make_k11(para[1:3])
  k11_r <- make_k11(para[4:6])
  f <- -0.5 * t(z) %*% ginv(k11) %*% z - 
    0.5 * log(det(k11)) +
    0.5 * m * log(2 * pi)
  f <- -f
  return(f)
}
#### Optimize
ln_likelihood(para)
a <- optimx(para, ln_likelihood, method = "BFGS")
## Calculate k21
make_k21 <- function(){
  k21_h <- make_cov(m1 = testing[c("x","y")], m2 = n_sws[c("x","y")],
                    kernel = kernel_fun, sof = sof_h, sd = sd)
  k21_v <- make_cov(m1 = testing["z"], m2 = n_sws["z"],
                    kernel = kernel_fun, sof = sof_v, sd = sd)
  k21 <- k21_h * k21_v
  return(k21)
}
k21 <- make_k21()

#Add noise
if (noise){
  r <- diag(nrow(k11))
  k11 <- k11 + r
}

## Predict
testing$nsws <- k21 %*% ginv(k11) %*% n_sws$nsws

## plot predicting
test_pic <- ggplot() +
  geom_point(data=n_sws,mapping = aes(x=nsws,y=z))+
  geom_line(data=testing, mapping = aes(x=nsws,y=z),color="#D3323F",orientation = "y")+
  scale_y_reverse()+
  facet_wrap(~x, nrow = 2)
test_pic
