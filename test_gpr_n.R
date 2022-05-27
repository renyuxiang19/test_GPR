library(tidyverse)
source("functions.r")

## Parameters of GPR
kernel_fun_name <- "Gaussian"
sof_h <- 5
sof_v <- 5
sd <- 20
nu <- 1.25

## read 2D data.
n_sws <- read_MAIC("kaminokoike_SWS.dat")
###  
n_sws$y <- 0

## Plot raw data.

raw_pic <- ggplot(data=n_sws) +
  geom_point(mapping = aes(x=nsws,y=z))+
  geom_line(mapping = aes(x=nsws,y=z),orientation = "y")+
  scale_y_reverse()+
  facet_wrap(~x, nrow = 2)

## 

