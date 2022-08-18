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
    ##
    testvar = NULL,
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
    create_mesh = function(size_v = 0.1, size_h = 0){
      # Prepare testing data (mesh).
      depth_raw <- self$n_sws |> 
        dplyr::group_by(x) |> 
        dplyr::summarise(min_depth = min(z), max_depth = max(z)) |>
        dplyr::ungroup() |>
        dplyr::arrange(x)
      if (size_h <= 0) {
        depth <- depth_raw
      }else{
        x_denser <- private$divide_KP(points = depth_raw$x, size = size_h)
        depth <- data.frame(x = x_denser,
                            min = private$lerp(depth_raw$x, depth_raw$min_depth, new_x = x_denser),
                            max = private$lerp(depth_raw$x, depth_raw$max_depth, new_x = x_denser))
        if (nrow(depth) > nrow(depth_raw)) {
          private$whether_contour <- TRUE
        }
      }
      self$testing <- as.list(depth) |> 
        rlang::set_names(NULL) |>
        purrr::pmap_dfr(function(dis, min, max){
          value <- seq(min, max, by = size_v)
          len <- length(value)
          tibble::tibble(x = rep(dis,times = len), z = value)
        }) |>
        mutate(y=0)
      private$whether_mesh <- TRUE
      cat("GPR: Mesh has been created.", "\n")
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
        names(self$para) <- private$para_name_nonu
        private$whether_nu <- FALSE
      }
      private$whether_set_parameter <- TRUE
      private$para_ini <- self$para
      cat("GPR: Parameters have been set.", "\n")
      invisible(self)
    },
    #
    set_par_scope = function(sof_h_t,sof_v_t, sd_t, sof_h_r, sof_v_r, sd_r, nu){
      # Check if nu is used.
      if (!rlang::is_missing(nu)) {
        self$para <- c(NA, NA, NA, NA, NA, NA, NA) 
        names(self$para) <- private$para_name
        par_scope <- list(sof_h_t, sof_v_t, sd_t, sof_h_r, sof_v_r, sd_r, nu)
        private$whether_nu <- TRUE
        private$kernel_fun <- "kernel_wm"
      }else{
        self$para <- c(NA, NA, NA, NA, NA, NA)
        names(self$para) <- private$para_name_nonu
        par_scope <- list(sof_h_t, sof_v_t, sd_t, sof_h_r, sof_v_r, sd_r)
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
      cat("GPR: The scope of paramters have been set and parameters (self$para) have been reset to be NA.", "\n")
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
      self$testing <- dplyr::mutate(self$testing,
                                    nsws = {k21 %*% ginv(k11) %*% self$n_sws$nsws} |> as.vector())
      # plot
      private$plot_predict()
      invisible(self)
    },
    opt = function(mode){
      stopifnot(whether_mesh)
      private$prepare_likelihood()
      match.arg(mode, c("BFGS", "GA"))
      switch (mode,
              BFGS = private$opt_bfgs(par = self$para, func = private$ln_likelihood),
              GA = private$opt_ga(lower = private$lower, upper = private$upper, func = private$ln_likelihood)
      )
      cat("Parameters have been optimized.", "\n")
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
    noise = TRUE,
    upper = NA,
    lower = NA,
    para_name = c("sof_h_t", "sof_v_t", "sd_t", "sof_h_r", "sof_v_r", "sd_r","nu"),
    para_name_nonu = c("sof_h_t", "sof_v_t", "sd_t", "sof_h_r", "sof_v_r", "sd_r"),
    para_ini = NULL,
    # Parameters used for likelihood function.
    z = NA,
    thirdterm = NA,
    #
    whether_mesh = FALSE,
    whether_nu = FALSE,
    whether_set_parameter = FALSE,
    whether_set_scope = FALSE,
    whether_contour = FALSE,
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
    plot_predict = function(...){
      if (private$whether_contour) {
        private$plot_predict_contour(data = self$testing, ...)
      }else{
        private$plot_predict_line(...)
      }
    },
    plot_predict_line = function(row = 2){
      self$test_pic <- ggplot() +
        geom_point(data = self$n_sws, mapping = aes(x = nsws, y = z))+
        geom_line(data = self$testing, mapping = aes(x = nsws, y = z),color="#D3323F",orientation = "y")+
        scale_y_reverse()+
        facet_wrap(~x, nrow = row)
      self$test_pic |> print()
    },
    plot_predict_contour = function(data, TitleColorbar, TitleX, TitleY, begin_X, max_Z,
                                    font_family = 'Times New Roman',font_size = 30, width = 1200 , height = 600){
      XYZtoZ <- function(rawdata){
        # this function is used to read and tidy data for contour plots (https://plotly.com/r/contour-plots/).
        # The input data should be a data frame with 3 columns, 
        # 1st and 2nd columns are for X and Y coordinates, respectively.
        # 3rd column stands for the Z value you want to plot.
        # this function convert the Z date to a 2D matrix for function plot_ly in the library(Plotly)
        # the the file of input data is supposed to be *.csv without header
        
        # Arguments
        # raw data : the input data frame
        
        colnames(rawdata)<-c("X","Y","Z")
        rawdata$X
        a<-tapply(rawdata$Y, rawdata$X, length)
        if(any(a!=a[1])){
          stop("Input data is not a full matrix")
        }
        xlen<-length(a)
        ylen<-a[1]
        zdata<-matrix(0, nrow = ylen, ncol = xlen)
        x<-as.numeric(attr(factor(rawdata$X),"levels"))
        y<-as.numeric(attr(factor(rawdata$Y),"levels"))
        zdalist<-split(rawdata$Z,rawdata$Y)
        for (i in 1:ylen) {
          zdata[i,]<-zdalist[[i]]
        }
        outda<-list(Z=zdata,X=x,Y=y)
        return(outda)
      }
      # tidy data
      DAT <- XYZtoZ(data)
      #
      fig_arguments <- list(
        title_colorbar = "N_sws",
        title_X = "Horizontal distance (m)",
        title_Y = "Depth (m)",
        shifX = 0
      )
      if(!missing(TitleColorbar)) fig_arguments$title_colorbar <- TitleColorbar
      if(!missing(TitleX)) fig_arguments$title_X <- TitleX
      if(!missing(TitleY)) fig_arguments$title_Y <- TitleY
      if(!missing(begin_X)) fig_arguments$shifX <- begin_X 
      # make the Y Axes to begin from 0 at the top
      Y_min<-min(DAT$Y)
      DAT$Y<- DAT$Y-Y_min
      DAT$Y<-rev(DAT$Y)   #reverse the value of Y axes
      Y_max<-max(DAT$Y)
      # shift the X axes
      shif_X <- fig_arguments$shifX
      DAT$X<-DAT$X+shif_X
      X_min<-min(DAT$X)
      # range of Z axes
      if(missing(max_Z)) max_Z <- max(DAT$Z)
      Z_max <- max_Z
      #
      fig <- plot_ly(
        type = 'contour',
        z = DAT$Z,
        y = DAT$Y,
        x = DAT$X,
        autocontour = F,
        colorscale= list(c(0, 0.33,0.66, 1), c('#FF0000', '#FFFB00', '#339502', '#023AB2')),
        contours = list(
          start = 0,
          end = Z_max,
          size = Z_max/10,
          showlines=FALSE
        ),
        colorbar=list(
          title=list(text=fig_arguments$title_colorbar,font=list(size=font_size)),
          tickfont=list(size=font_size),
          len=0.9,
          exponentformat="power"    
        ),
        line = list(smoothing = 0)  
      )
      font_for_tick <- list(        
        size = font_size
      )
      x <- list(
        title = fig_arguments$title_X,          
        titlefont = list(size=font_size),  
        autotick = FALSE,
        ticks = "outside",
        tickmode = "linear",
        tick0 = X_min,            
        dtick = 5,                 
        ticklen = 8,              
        tickfont = font_for_tick,   
        tickwidth = 1,              
        tickcolor = toRGB("black")
      )
      y <- list(
        title=fig_arguments$title_Y,
        titlefont = list(size=font_size),  
        ticklen = 5,
        tickfont = font_for_tick,  
        tickwidth = 1,
        range = c(Y_max,0)   
      )
      self$test_pic <- fig %>% layout(xaxis = x,yaxis = y,
                                      font = list(family = font_family ,size= font_size))
    },
    # ------ optimizing functions ------
    prepare_likelihood = function(){
      private$z <- self$n_sws$nsws |> matrix()
      m <- length(private$para)
      private$thirdterm <- 0.5 * m * log(2 * pi)
    },
    ln_likelihood = function(para){
      make_k <- private$createfunc_make_k(private$whether_nu)
      kernel_function <- get(private$kernel_fun, envir = private)
      k11_t <- make_k(para[-c(4:6)], cm1 = self$n_sws, cm2 = self$n_sws, kernel = kernel_function)
      k11_r <- make_k(para[-c(1:3)], cm1 = self$n_sws, cm2 = self$n_sws, kernel = kernel_function)
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
    #
    opt_bfgs = function(par, func){
      if (!private$whether_set_parameter) {
        stop("Error in GPR: pleas set initial parameters before using BFGS to optimize them. (using 'set_parameter')")
      }
      func <- private$more_info(func = func)
      opt <- optimx(par, func, method = "BFGS")
    },
    #
    opt_ga = function(lower, upper, func){
      if (!private$whether_set_scope) {
        stop("Error in GPR: pleas set the scope of each parameter before using GA to optimize them. (using 'set_par_scope')")
      }
      out_ga2 <- GA::ga(type = "real-valued", fitness = func, 
                        lower = lower, upper = upper,
                        popSize = 120, maxiter = 100, run = 20, parallel = 7,
                        optim = FALSE)
      name_para <- names(self$para)
      self$para <- out_ga2@solution[1,] |> as.numeric()
      names(self$para) <- name_para
      private$whether_set_parameter <- TRUE
    },
    #
    more_info = function(func){
      # used for likelihood function. 
      # Write the argument into 'self$para' and record the time the function has been run.
      force(func)
      time <- 0
      function(x, ...){
        res <- func(x, ...)
        names(res) <- NULL
        self$para <- x
        time <<- time + 1
        cat("BFGS | iter = ", time, " | The parameters are: ", "\n" ,sep = "")
        print(x)
        cat("\n")
        res
      }
    }
  )
)
