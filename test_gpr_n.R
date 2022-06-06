# Test GPR with real N-value data.
# library(tidyverse)
# library(rlang)
library(MASS)
library(optimx)
source("functions.r")

## (initial) Parameters of GPR.
filename <- "kaminokoike_SWS.dat"
kernel_fun <- "kernel_g"
sof_h_t <- 200   
sof_v_t <- 2 
sd_t <- 10
sof_h_r <- 0.01
sof_v_r <- 0.01
sd_r <- 10 
nu <- 1.25
noise <- TRUE
### Parameter vector.
para <- c(sof_h_t, sof_v_t, sd_t, sof_h_r, sof_v_r, sd_r)

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
print(raw_pic)
writeLines("Plot the raw data.")

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

#### Create cov matrix K11 or K21.
make_k <- function(par, cm1, cm2){
  # "cm1" and "cm2" are coordinate matrices with "x", "y", "z" as their column name.
  #   "x" and "y" are horizontal coordinate, "z" is depth.
  # "par" is a numerical vector for kernel function.
  sof_h <- par[1]
  sof_v <- par[2]
  sd <- par[3]
  #
  k21_h <- make_cov(m1 = cm1[c("x","y")], m2 = cm2[c("x","y")],
                    kernel = kernel_fun, sof = sof_h, sd = sd)
  k21_v <- make_cov(m1 = cm1["z"], m2 = cm2["z"],
                    kernel = kernel_fun, sof = sof_v, sd = sd)
  k21 <- k21_h * k21_v
  return(k21)
}

#### likelihood function.
ln_likelihood <- function(para){
  # Observation vector.
  z <- n_sws$nsws |> matrix()
  #
  m <- length(para)
  k11_t <- make_k(para[1:3], cm1 = n_sws, cm2 = n_sws)
  k11_r <- make_k(para[4:6], cm1 = n_sws, cm2 = n_sws)
  k11_t <- k11_t/k11_t[1,1]
  k11_r <- k11_r/k11_r[1,1]
  k11 <- k11_t + k11_r
  f <- -0.5 * t(z) %*% ginv(k11) %*% z - 
    0.5 * log(det(k11)) +
    0.5 * m * log(2 * pi)
  f <- -f
  f <- as.numeric(f)
  return(f)
}

#### Optimize parameters
opt_para <- function(par, func){
  func <- match.fun(func)
  len_par <- length(par)
  opt <- optimx(para, func, method = "BFGS")
  new_para <- opt[seq_len(len_par)] |> as.numeric()
  return(new_para)
}

writeLines("Optimize parameters...")
para <- opt_para(para, "ln_likelihood")

#### Calculate covariance matrices
k11 <- make_k(para[1:3], cm1 = n_sws, cm2 = n_sws) 
k21 <- make_k(para[1:3], cm1 = testing, cm2 = n_sws) 
#### Add noise to K11
if (noise){
  r <- diag(nrow(k11))
  k11 <- k11 + r
}

## Predict
testing$nsws <- k21 %*% ginv(k11) %*% n_sws$nsws

## plot predicting
test_pic <- ggplot() +
  geom_point(data = n_sws, mapping = aes(x = nsws, y = z))+
  geom_line(data = testing, mapping = aes(x = nsws, y = z),color="#D3323F",orientation = "y")+
  scale_y_reverse()+
  facet_wrap(~x, nrow = 2)
test_pic |> print()
writeLines("Plot predicting curve.")

