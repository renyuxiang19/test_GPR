library(MASS)
library(optimx)
library(GA)
library(doParallel)
library(foreach)
library(iterators)
library(R6)

GPR <- R6::R6Class(
  classname = "GPR",
  public = list(
    #variables
    para = NULL,
    ## Data
    n_sws = NULL,
    testing = NULL,
    ## Fig objects
    raw_pic = NULL,
    test_pic = NULL,
    #
    initialize = function(data_file){
      if(rlang::is_missing(data_file)){
        stop("Error: GPR: please specify a file (MAIC) to read data.")
      }
      # Read data.
      if (is.character(data_file)){
        self$n_sws = private$read_MAIC(data_file) |> 
          dplyr::mutate(y=0) |>
          tibble::rowid_to_column("ID") |>
          dplyr::select(-ID)|>
          dplyr::filter(-1 < x & x < 200)
        # Plot raw data.
        self$raw_pic = ggplot(data=self$n_sws)+
          geom_point(mapping = aes(x=nsws,y=z))+
          geom_line(mapping = aes(x=nsws,y=z),orientation = "y")+
          scale_y_reverse()+
          facet_wrap(~x, nrow = 2)
        print(self$raw_pic)
        writeLines("GPR: Plot the raw data.")
      }else{stop("Error: GPR: please use the name (string) to specify a file (MAIC).")}
    },
    #
    create_mesh = function(){
      # Prepare testing data (mesh).
      private$depth <- self$n_sws |> 
        group_by(x) |> 
        summarise(min_depth = min(z), max_depth = max(z))
      self$testing <- list(private$depth$x, private$depth$min_depth, private$depth$max_depth) |> 
        purrr::pmap_dfr(function(dis, min, max){
          value <- seq(min, max, by = private$mesh_size_v)
          len <- length(value)
          tibble::tibble(x = rep(dis,times = len), z = value)
        }) |>
        mutate(y=0)
    },
    #
    set_kernel = function(fun_name){
      if (is.character(fun_name)){
        private$kernel_fun <- fun_name
      }else{
        stop("Error: GPR: The name of kernel function should be a character string")
      }
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
        private$nu <- nu
        self$para <- c(sof_h_t, sof_v_t, sd_t,sof_h_r, sof_v_r, sd_r, nu) |> matrix(nrow = 1)
        colnames(self$para) <- para_name
        private$whether_nu <- TRUE
        private$kernel_fun <- "kernel_wm"
      }else{
        self$para <- c(sof_h_t, sof_v_t, sd_t,sof_h_r, sof_v_r, sd_r)
        colnames(self$para) <- c("sof_h_t", "sof_v_t", "sd_t", "sof_h_r", "sof_v_r", "sd_r")
        private$whether_nu <- FALSE
      }
      private$whether_set_parameter <- TRUE
    },
    #
    set_par_scope = function(sof_h_t,sof_v_t, sd_t, sof_h_r, sof_v_r, sd_r, nu){
      # Check if nu is used.
      if (!rlang::is_missing(nu)) {
        par_scope <- list(sof_h_t, sof_v_t, sd_t,sof_h_r, sof_v_r, sd_r, nu)
        private$whether_nu <- TRUE
        private$kernel_fun <- "kernel_wm"
      }else{
        par_scope <- list(sof_h_t, sof_v_t, sd_t,sof_h_r, sof_v_r, sd_r)
        private$whether_nu <- FALSE
      }
      # Check the argument format.
      if (!purrr::every(par_scope, ~length(.x)== 2)) {
        stop("Error: GPR: The scope of parameters should be vectors which include upper and lower limit. like nu = c(0, 1)")
      }
      #
      private$upper <- par_scope |> purrr::map_dbl(max)
      private$lower <- par_scope |> purrr::map_dbl(min)
      private$whether_set_scope <- TRUE
    },
    
    
  ),
  #
  private = list(
    kernel_fun = "kernel_g",
    sof_h_t = NULL,
    sof_v_t = NULL,
    sd_t = NULL,
    sof_h_r = NULL,
    sof_v_r = NULL,
    sd_r = NULL ,
    nu = NULL,
    depth = NULL,
    noise = TRUE,
    mesh_size_v = 0.1,
    upper = NULL,
    lower = NULL,
    para_name = c("sof_h_t", "sof_v_t", "sd_t", "sof_h_r", "sof_v_r", "sd_r","nu"),
    # Parameters used for likelihood function.
    z = NULL,
    m = NULL,
    thirdterm = NULL,
    #
    whether_nu = FALSE,
    whether_set_parameter = FALSE,
    whether_set_scope = FALSE,
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
      ############################################################
      kernel <- match.fun(kernel)
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
    #
    make_k = function(par, cm1, cm2){
      # "cm1" and "cm2" are coordinate matrices with "x", "y", "z" as their column name.
      #   "x" and "y" are horizontal coordinate, "z" is depth.
      # "par" is a numerical vector for kernel function.
      sof_h <- par[1]
      sof_v <- par[2]
      sd <- par[3]
      nu <- par[4]
      #
      k21_h <- private$make_cov(m1 = cm1[c("x","y")], m2 = cm2[c("x","y")],
                        kernel = self$kernel_fun, sof = sof_h, 
                        sd = sd, nu = nu)
      k21_v <- private$make_cov(m1 = cm1["z"], m2 = cm2["z"],
                        kernel = self$kernel_fun, sof = sof_v, 
                        sd = sd, nu = nu)
      k21 <- k21_h * k21_v
      return(k21)
    },
    prepare_likelihood = function(){
      private$z <- self$n_sws$nsws |> matrix()
      private$m <- length(private$para)
      private$thirdterm <- 0.5 * private$m * log(2 * pi)
    }
  )
)
