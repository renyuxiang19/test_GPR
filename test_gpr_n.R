# Test GPR with real N-value data.
library(tidyverse)
library(rlang)
# library(GenSA)
library(MASS)
library(optimx)
library(GA)
library(doParallel)
library(foreach)
library(iterators)
source("functions.r")

## (initial) Parameters of GPR.
filename <- "kaminokoike_SWS.dat"
kernel_fun <- "kernel_wm"
sof_h_t <- 30
sof_v_t <- 5
sd_t <- 2
sof_h_r <- 0.01
sof_v_r <- 0.01
sd_r <- 5
nu <- 5
noise <- TRUE
### Parameter vector.
para <- c(sof_h_t, sof_v_t, sd_t,sof_h_r, sof_v_r, sd_r, nu)
lower <- c( 1,0.1, 1, 1,0.1, 1, 0.5)
upper <- c(20, 10,50,20, 10,50, 10)
### Unimportant parameter
mesh_size_v <- 0.25
mesh_size_h <- 10

## read 2D data.
n_sws <- read_MAIC(filename) |> 
  dplyr::mutate(y=0) |>
  tibble::rowid_to_column("ID") |>
  dplyr::select(-ID)|>
  dplyr::filter(-1 < x & x < 200)

## Plot raw data.
raw_pic <- ggplot(data=n_sws)+
  geom_point(mapping = aes(x=nsws,y=z))+
  geom_line(mapping = aes(x=nsws,y=z),orientation = "y")+
  scale_y_reverse()+
  facet_wrap(~x, nrow = 2)
print(raw_pic)
writeLines("Plot the raw data.")

## Prepare testing data (mesh).
depth_raw <- n_sws |> 
  dplyr::group_by(x) |> 
  dplyr::summarise(min_depth = min(z), max_depth = max(z)) |>
  dplyr::ungroup() |>
  dplyr::arrange(x)
x_denser <- divide_KP(points = depth_raw$x, size = mesh_size_h)
depth <- data.frame(x = x_denser,
                    min = lerp(depth_raw$x, depth_raw$min_depth, new_x = x_denser),
                    max = lerp(depth_raw$x, depth_raw$max_depth, new_x = x_denser))
testing <- as.list(depth) |> 
  rlang::set_names(NULL) |>
  purrr::pmap_dfr(function(dis, min, max){
    value <- seq(min, max, by = mesh_size_v)
    len <- length(value)
    tibble::tibble(x = rep(dis,times = len), z = value)
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
  nu <- par[4]
  #
  k21_h <- make_cov(m1 = cm1[c("x","y")], m2 = cm2[c("x","y")],
                    kernel = kernel_fun, sof = sof_h, sd = sd, nu = nu)
  k21_v <- make_cov(m1 = cm1["z"], m2 = cm2["z"],
                    kernel = kernel_fun, sof = sof_v, sd = sd, nu = nu)
  k21 <- k21_h * k21_v
  return(k21)
}
#### likelihood function.
  # Observation vector.
z <- n_sws$nsws |> matrix()
m <- length(para)
thirdterm <- 0.5 * m * log(2 * pi)
#
ln_likelihood <- function(para){
  k11_t <- make_k(para[-c(4:6)], cm1 = n_sws, cm2 = n_sws)
  k11_r <- make_k(para[-c(1:3)], cm1 = n_sws, cm2 = n_sws)
  k11_r_pivot <- k11_r[1,1]
  k11_t <- k11_t/k11_r_pivot
  k11_r <- k11_r/k11_r_pivot
  k11 <- k11_t + k11_r
  k11 <- `diag<-`(k11, diag(k11) + k11[1,1]*0.1)
  f <- -0.5 * t(z) %*% ginv(k11*k11_r_pivot) %*% z - 
    0.5 *(log(k11_r_pivot) * nrow(k11)+ log(det(k11))) + thirdterm
  #f <- -f
  f <- as.numeric(f)
  return(f)
}
#
ln_likelihood_2 <- function(para){
  k11_t <- make_k(para, cm1 = n_sws, cm2 = n_sws)
  k11_t_pivot <- k11_t[1,1]
  k11_t <- k11_t/k11_t_pivot
  k11_t <- `diag<-`(k11_t, diag(k11_t) + 0.1)
  k11 <- k11_t
  f <- -0.5 * t(z) %*% ginv(k11*k11_t_pivot) %*% z - 
    0.5 *(log(k11_t_pivot) * nrow(k11) + log(det(k11))) + thirdterm
  #f <- -f
  f <- as.numeric(f)
  return(f)
}

#### Optimize parameters
writeLines("Optimize parameters...")
opt_para <- function(par, func){
  func <- match.fun(func)
  len_par <- length(par)
  opt <- optimx(para, func, method = "BFGS")
  new_para <- opt[seq_len(len_par)] |> as.numeric()
  return(new_para)
}
# out_sa <- GenSA::GenSA(par = para, fn = ln_likelihood, 
#              lower = lower, upper = upper,
#              control = list(smooth = TRUE , max.call = 14))
out_ga2 <- GA::ga(type = "real-valued", fitness = ln_likelihood, 
                 lower = lower, upper = upper,
                 popSize = 120, maxiter = 100, run = 20, parallel = 7,
                 optim = FALSE)
para <- out_ga2@solution[1,] |> as.numeric()
para_out <- out_ga2@solution
colnames(para_out) <- c("sof_h_t", "sof_v_t", "sd_t", "sof_h_r", "sof_v_r", "sd_r","nu")
# out <- optimization::optim_sa(fun = ln_likelihood, start = para,
#                               lower = lower, upper = upper,
#                               control = list(nlimit = 2))
#para <- opt_para(para, "ln_likelihood")

#### Calculate covariance matrices
k11 <- make_k(para[-c(4:6)], cm1 = n_sws, cm2 = n_sws) 
# k11_r <- make_k(para[-c(1:3)], cm1 = n_sws, cm2 = n_sws) 
k21 <- make_k(para[-c(4:6)], cm1 = testing, cm2 = n_sws) 
# k21_r <- make_k(para[-c(1:3)], cm1 = testing, cm2 = n_sws) 
k11 <- k11 #+ k11_r
k21 <- k21 #+ k21_r
#### Add noise to K11
if (noise){
  k11 <- `diag<-`(k11, diag(k11)+1)
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
#
writeLines("Optimized Parameters:")
print(para_out)

