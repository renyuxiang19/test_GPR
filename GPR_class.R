library(MASS)
library(optimx)
library(GA)
library(doParallel)
library(foreach)
library(iterators)
library(R6)
library(tidyverse)

GPR <- R6::R6Class(
  classname = "GPR",
  public = list(
    #variables
    para = NA,
    ## Data
    n_sws = NA,
    testing = NA,
    ## Fig objects
    raw_pic = NA,
    test_pic = NA,
    #
    initialize = function(data_file){
      if(rlang::is_missing(data_file)){
        stop("Error in GPR: please specify a file (MAIC) to read data.")
      }
      # Read data.
      if (is.character(data_file)){
        self$n_sws = private$read_MAIC(data_file) |> 
          dplyr::mutate(y=0) |>
          tibble::rowid_to_column("ID") |>
          dplyr::select(-ID)|>
          dplyr::filter(-1 < x & x < 200)
        # Plot raw data.
        self$raw_pic = ggplot2::ggplot(data=self$n_sws)+
          geom_point(mapping = aes(x=nsws,y=z))+
          geom_line(mapping = aes(x=nsws,y=z),orientation = "y")+
          scale_y_reverse()+
          facet_wrap(~x, nrow = 2)
        print(self$raw_pic)
        writeLines("GPR: Plot the raw data.")
      }else{stop("Error in GPR: please use the name (string) to specify a file (MAIC).")}
      invisible(self)
    },
    #
    create_mesh = function(size_v = 0.1){
      # Prepare testing data (mesh).
      private$depth <- self$n_sws |> 
        group_by(x) |> 
        summarise(min_depth = min(z), max_depth = max(z))
      self$testing <- list(private$depth$x, private$depth$min_depth, private$depth$max_depth) |> 
        purrr::pmap_dfr(function(dis, min, max){
          value <- seq(min, max, by = size_v)
          len <- length(value)
          tibble::tibble(x = rep(dis,times = len), z = value)
        }) |>
        mutate(y=0)
      private$whether_mesh <- TRUE
      invisible(self)
    },
    #
    set_kernel = function(fun_name){
      if (is.character(fun_name)){
        private$kernel_fun <- fun_name
      }else{
        stop("Error in GPR: The name of kernel function should be a character string")
      }
      invisible(self)
    },
    #
    set_parameter = function(sof_h_t,sof_v_t, sd_t, sof_h_r, sof_v_r, sd_r, nu){
      private$sof_h_t <- sof_h_t
      private$sof_v_t <- sof_v_t
      private$sd_t <- sd_t
      private$sof_h_r <- sof_h_r
      private$sof_v_r <- sof_v_r
      private$sd_r <- sd_r
      if (!rlang::is_missing(nu)) {
        self$para <- c(sof_h_t, sof_v_t, sd_t,sof_h_r, sof_v_r, sd_r, nu) 
        names(self$para) <- private$para_name
        private$whether_nu <- TRUE
        private$kernel_fun <- "kernel_wm"
      }else{
        self$para <- c(sof_h_t, sof_v_t, sd_t,sof_h_r, sof_v_r, sd_r)
        names(self$para) <- c("sof_h_t", "sof_v_t", "sd_t", "sof_h_r", "sof_v_r", "sd_r")
        private$whether_nu <- FALSE
      }
      private$whether_set_parameter <- TRUE
      private$para_ini <- self$para
      invisible(self)
    },
    #
    set_par_scope = function(sof_h_t,sof_v_t, sd_t, sof_h_r, sof_v_r, sd_r, nu){
      # Check if nu is used.
      if (!rlang::is_missing(nu)) {
        self$para <- c(NA, NA, NA,NA, NA, NA, NA) 
        names(self$para) <- private$para_name
        par_scope <- list(sof_h_t, sof_v_t, sd_t,sof_h_r, sof_v_r, sd_r, nu)
        private$whether_nu <- TRUE
        private$kernel_fun <- "kernel_wm"
      }else{
        self$para <- c(NA, NA, NA,NA, NA, NA, NA)
        names(self$para) <- c("sof_h_t", "sof_v_t", "sd_t", "sof_h_r", "sof_v_r", "sd_r")
        par_scope <- list(sof_h_t, sof_v_t, sd_t,sof_h_r, sof_v_r, sd_r)
        private$whether_nu <- FALSE
      }
      # Check the argument format.
      if (!purrr::every(par_scope, ~length(.x)== 2)) {
        stop("Error in GPR: The scope of parameters should be vectors which include upper and lower limit. like nu = c(0, 1)")
      }
      #
      private$upper <- par_scope |> purrr::map_dbl(max)
      private$lower <- par_scope |> purrr::map_dbl(min)
      private$whether_set_scope <- TRUE
      invisible(self)
    },
    predict = function(){
      private$check_predict()
      make_k <- private$createfunc_make_k(private$whether_nu)
      kernel_function <- get(private$kernel_fun, envir = private)
      # Calculate covariance matrices of trend component.
      k11 <- make_k(self$para[-c(4:6)], cm1 = self$n_sws, cm2 = self$n_sws, kernel = kernel_function) 
      k21 <- make_k(self$para[-c(4:6)], cm1 = self$testing, cm2 = self$n_sws, kernel = kernel_function) 
      if (private$noise){
        k11 <- `diag<-`(k11, diag(k11) + 1)
      }
      # predict
      self$testing$nsws <- k21 %*% ginv(k11) %*% self$n_sws$nsws
      # plot
      private$plot_predict()
      invisible(self)
    }
  ),
  #
  private = list(
    kernel_fun = "kernel_g",
    sof_h_t = NA,
    sof_v_t = NA,
    sd_t = NA,
    sof_h_r = NA,
    sof_v_r = NA,
    sd_r = NA ,
    nu = NA,
    depth = NA,
    noise = TRUE,
    upper = NA,
    lower = NA,
    para_name = c("sof_h_t", "sof_v_t", "sd_t", "sof_h_r", "sof_v_r", "sd_r","nu"),
    para_ini = NULL,
    # Parameters used for likelihood function.
    z = NA,
    thirdterm = NA,
    #
    whether_mesh = FALSE,
    whether_nu = FALSE,
    whether_set_parameter = FALSE,
    whether_set_scope = FALSE,
    ######################### Pure Functions ##############################
    # Read data from a file. 
    read_MAIC = function(file_name){
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
        if (!private$is_mathinteger(data_number[i])) {
          paste("Error in read_MAIC: the number of data at ", i,
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
    },
    # Kernel functions
    kernel_g = function(d,sof,sd){
      k <- (sd^2) * exp(-pi* (d/sof)^2) 
      return(k)
    },
    kernel_m = function(d,sof,sd){
      k <- (sd^2) * exp(-2* (d/sof))
      return(k)
    },
    kernel_b = function(d,sof,sd){
      dif <- d
      if (dif < sof) {
        k <- (sd^2) * (1-(dif/sof)) 
      }else{
        k <- 0
      }
      return(k)
    },
    kernel_wm = function(d, nu, sof, sd){
      if (d == 0) {
        covariance <- sd^2
      }else{
        frac <- sqrt(pi)*gamma(nu+0.5)*d/(gamma(nu)*sof)
        covariance <- sd^2 * (2/gamma(nu)) * frac^nu * besselK(2*frac, nu)
      }
      return(covariance)
    },
    # Calculate covariance matrix
    make_cov = function(m1, m2, kernel, ... ){
      # "m1" and "m2" should be matrices or data.frames with same numbers of columns.
      #   e.g., if you want to calculate in two horizontal directions, the "m1" and "m1" should both have 2 columns.
      #   The first column in the both two matrices means one direction, and the second for another.
      #
      #kernel <- match.fun(kernel)
      # Check input data.
      if (ncol(m1) != ncol(m2)){ stop("Error: The number of columns of input matrices are not consistent.") }
      # Create distance matrices
      d <- rdist::cdist(m1, m2) 
      # Calculate covariance matrices with kernel function
      k <- purrr::modify(d, kernel, ...)
      return(k)
    },
    #
    is_mathinteger = function(v){
      all(round(v) == v)
    },
    createfunc_make_k = function(use_nu){
      # "cm1" and "cm2" are coordinate matrices with "x", "y", "z" as their column name.
      #   "x" and "y" are horizontal coordinate, "z" is depth.
      # "par" is a numerical vector for kernel function.
      if(use_nu){
        function(par, cm1, cm2, kernel){
          sof_h <- par[1]
          sof_v <- par[2]
          sd <- par[3]
          nu <- par[4]
          #
          if(is.character(kernel)){
            kernel <- get(kernel, envir = private)
          }
          #
          k21_h <- private$make_cov(m1 = cm1[c("x","y")], m2 = cm2[c("x","y")],
                                    kernel = kernel, sof = sof_h, 
                                    sd = sd, nu = nu)
          k21_v <- private$make_cov(m1 = cm1["z"], m2 = cm2["z"],
                                    kernel = kernel, sof = sof_v, 
                                    sd = sd, nu = nu)
          k21 <- k21_h * k21_v
          return(k21)
        }
      }else{
        function(par, cm1, cm2, kernel){
          sof_h <- par[1]
          sof_v <- par[2]
          sd <- par[3]
          #
          if(is.character(kernel)){
            kernel <- get(kernel, envir = private)
          }
          #
          k21_h <- private$make_cov(m1 = cm1[c("x","y")], m2 = cm2[c("x","y")],
                                    kernel = kernel, sof = sof_h, 
                                    sd = sd)
          k21_v <- private$make_cov(m1 = cm1["z"], m2 = cm2["z"],
                                    kernel = kernel, sof = sof_v, 
                                    sd = sd)
          k21 <- k21_h * k21_v
          return(k21)
        }
      }
    },
    #
    divide_KP = function(points, size){
      # divide by key points. Those key points will be kept.
      # find the intervals.
      points_lag <- points[-length(points)]
      points <- points[-1]
      interval <- {points - points_lag} |> abs()
      # divide each interval and round it up.
      division <- {interval / size} |> ceiling()
      
      # calculate points in each interval.
      points_new <- list(points_lag, points, division) |>
        purrr::pmap(function(start, end, number){
          seq(start, end, length = number + 1)
        }) |> unlist() |> unique()
      return(points_new)
    },
    #
    lerp = function(x, y, new_x){
      # Linear interpolation
      interp_f <- approxfun(x, y)
      new_y <- interp_f(new_x)
      return(new_y)
    },
    ####################### End of pure functions ##############################
    prepare_likelihood = function(){
      private$z <- self$n_sws$nsws |> matrix()
      m <- length(private$para)
      private$thirdterm <- 0.5 * m * log(2 * pi)
    },
    #
    check_predict = function(){
     if(purrr::some(self$para, is.na)){
       stop("Error in GPR: pleas use 'set_parameter' method to set ALL parameters (set nu only when you want to use WM model).")
     }
     if(private$kernel_fun == "kernel_wm"){
       if (private$whether_nu  == FALSE) {
         stop("Error in GPR: you have not set the initial value of parameter nu.")
       }
     }else{
       if (private$whether_nu  == TRUE) {
         warning("GPR: you set a value of nu without using WM model.")
       }
     }
     if(!private$whether_mesh){
       stop("Error in GPR: please use 'create_mesh' to creat mesh.")
     }
    },
    match_kernel = function(name){
    },
    # ------ plotting functions ------
    #
    plot_predict = function(){
      self$test_pic <- ggplot() +
        geom_point(data = self$n_sws, mapping = aes(x = nsws, y = z))+
        geom_line(data = self$testing, mapping = aes(x = nsws, y = z),color="#D3323F",orientation = "y")+
        scale_y_reverse()+
        facet_wrap(~x, nrow = 2)
      self$test_pic |> print()
    },
    # ------ optimizing functions ------
    ln_likelihood = function(para){
      make_k <- createfunc_make_k(private$whether_nu)
      k11_t <- make_k(para[-c(4:6)], cm1 = self$n_sws, cm2 = self$n_sws)
      k11_r <- make_k(para[-c(1:3)], cm1 = self$n_sws, cm2 = self$n_sws)
      k11_r_pivot <- k11_r[1,1]
      k11_t <- k11_t/k11_r_pivot
      k11_r <- k11_r/k11_r_pivot
      k11 <- k11_t + k11_r
      k11 <- `diag<-`(k11, diag(k11) + k11[1,1]*0.1)
      f <- -0.5 * t(private$z) %*% ginv(k11*k11_r_pivot) %*% private$z - 
        0.5 *(log(k11_r_pivot) * nrow(k11)+ log(det(k11))) + private$thirdterm
      #f <- -f
      f <- as.numeric(f)
      return(f)
    },
    opt_bfgs = function(par, func){
      len_par <- length(par)
      opt <- optimx(para, func, method = "BFGS")
      new_para <- opt[seq_len(len_par)] |> as.numeric()
      return(new_para)
    },
    opt_ga = function(lower, upper, func){
      out_ga2 <- GA::ga(type = "real-valued", fitness = func, 
                        lower = lower, upper = upper,
                        popSize = 120, maxiter = 100, run = 20, parallel = 7,
                        optim = FALSE)
      self$para <- out_ga2@solution[1,] |> as.numeric()
    }
    opt = function(mode){
      
    }
  )
)
