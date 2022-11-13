library(GA)
library(doParallel)
library(foreach)
library(iterators)
library(R6)
library(tidyverse)
library(plotly)
library(interp)
library(MASS)
#
# How to use: 
# 'new("file name")' -> 'set_parameter' or 'set_par_scope'  -> 'create_mesh'
# -> 'opt' -> 'predict'
GPR <- R6::R6Class(
  classname = "GPR",
  public = list(
    #variables
    para = NA,
    ## Data
    n_sws = NA,
    testing = NA,
    normalized_n = NULL,
    ## Fig objects
    raw_pic = NA,
    test_pic = NA,
    norm_pic = NA,
    ## opt result
    result_ga = NULL,
    result_BFGS = NULL,
    ##
    #
    initialize = function(data_file, plot_silent = FALSE){
      if(rlang::is_missing(data_file)){
        stop("Error in GPR: please specify a file (MAIC) to read data.")
      }
      self$plot_silent(plot_silent)
      # Read data.
      if (is.character(data_file)){
        self$n_sws = private$read_MAIC(data_file) |> 
          dplyr::mutate(y=0) |>
          tibble::rowid_to_column("ID") |>
          dplyr::select(-ID)
        private$raw_n_sws <- self$n_sws
        # Plot raw data.
        private$plot_raw()
        writeLines("GPR: Data has been read. Please create a mesh.")
      }else{stop("Error in GPR: please use the name (string) to specify a file (MAIC).")}
      invisible(self)
    },
    #
    create_mesh = function(size_v = private$mesh_size_v, 
                           size_h = private$mesh_size_h, 
                           rigid = FALSE, from_zero = FALSE, silent = FALSE){
      if (rigid) {
        grDevices::dev.off()
        self$n_sws <- private$rigidify(dat = private$raw_n_sws, size_v = size_v, 
                                       from_zero = from_zero)
        private$plot_raw()
        private$whether_rigid <- TRUE
      }else{
        self$n_sws <- private$raw_n_sws
      }
      self$testing <- private$mesh_data(size_v = size_v, size_h = size_h)
      private$whether_mesh <- TRUE
      if(!silent){
        cat("GPR: Mesh has been created.", "\n")
        if (rigid) {
        cat("     Using rigid mode, the depth of data may be modified.", "\n")
        }
      }
      invisible(self)
    },
    set_mesh_size = function(size_v, size_h){
      private$mesh_size_v <- size_v
      private$mesh_size_h <- size_h
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
        cat("GPR: WM model will be used.", "\n")
      }else{
        self$para <- c(sof_h_t, sof_v_t, sd_t,sof_h_r, sof_v_r, sd_r)
        names(self$para) <- private$para_name_nonu
        private$whether_nu <- FALSE
        private$kernel_fun <- "kernel_g"
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
        if (!private$whether_set_parameter) { 
          self$para <- c(NA, NA, NA, NA, NA, NA, NA)
          names(self$para) <- private$para_name
        }
        par_scope <- list(sof_h_t, sof_v_t, sd_t, sof_h_r, sof_v_r, sd_r, nu)
        private$whether_nu <- TRUE
        private$kernel_fun <- "kernel_wm"
        cat("GPR: WM model will be used.", "\n")
      }else{
        if (!private$whether_set_parameter) {
          self$para <- c(NA, NA, NA, NA, NA, NA)
          names(self$para) <- private$para_name_nonu
        }
        par_scope <- list(sof_h_t, sof_v_t, sd_t, sof_h_r, sof_v_r, sd_r)
        private$whether_nu <- FALSE
        private$kernel_fun <- "kernel_g"
      }
      # Check the argument format.
      if (!purrr::every(par_scope, ~length(.x)== 2)) {
        stop("Error in GPR: The scope of parameters should be vectors which include upper and lower limit. like nu = c(0, 1)")
      }
      #
      private$upper <- par_scope |> purrr::map_dbl(max)
      private$lower <- par_scope |> purrr::map_dbl(min)
      private$whether_set_scope <- TRUE
      cat("GPR: The scope of parameters have been set.", "\n")
      invisible(self)
    },
    #
    set_limit_x = function(scope){
      if (is.numeric(scope) & length(scope) == 2) {
        private$limit_x[1] <- min(scope)
        private$limit_x[2] <- max(scope)
        self$n_sws <-  private$raw_n_sws |> 
          dplyr::filter(private$limit_x[1] < x & x < private$limit_x[2])
        private$plot_raw()
        #
        if (private$whether_mesh == TRUE) {
          private$whether_mesh <- FALSE
          cat("GPR: Limit of x has been set. Please create a new mesh once again.", "\n")
        }else{cat("GPR: Limit of x has been set. Please create a mesh.", "\n")}
      }else{
        cat("Error in GPR: please set the upper and lower limit of x by a length 2 numerical vector.", "\n")
      }
      invisible(self)
    },
    #
    predict = function(){
      private$check_predict() # Check the condition for prediction.
      make_k <- private$createfunc_make_k(private$whether_nu)   # This is a function to create K matrix.
      kernel_function <- get(private$kernel_fun, envir = private)  # get the kernel function according to its name.
      # Calculate covariance matrices of trend component.
      k11 <- make_k(self$para[-c(4:6)], cm1 = self$n_sws, cm2 = self$n_sws, kernel = kernel_function) 
      k21 <- make_k(self$para[-c(4:6)], cm1 = self$testing, cm2 = self$n_sws, kernel = kernel_function) 
      if (private$noise){                      #add noise
        k11 <- `diag<-`(k11, diag(k11) + 1)
      }
      # predict
      self$testing <- dplyr::mutate(self$testing,
                                    nsws = {k21 %*% ginv(k11) %*% self$n_sws$nsws} |> as.vector())
      # plot
      self$plot_predict()
      # output log on screen.
      private$repdict_log()
      #
      private$whether_predict <- TRUE
      invisible(self)
    },
    opt = function(mode = "GA"){
      stopifnot(private$whether_mesh)
      private$prepare_likelihood()
      match.arg(mode, c("BFGS", "GA"))
      cat("GPR: Start to optimize parameters by ",mode," ...",sep = "", "\n")
      switch (mode,
              BFGS = private$opt_bfgs(par = self$para, func = private$ln_likelihood,
                                      lower = private$lower, upper = private$upper),
              GA = private$opt_ga(lower = private$lower, upper = private$upper, func = private$ln_likelihood)
      )
      cat("GPR: Parameters have been optimized.", "\n")
      invisible(self)
    },
    set_plot_lines = function(L = TRUE){
      private$whether_contour <- !L
      invisible(self)
    },
    plot_predict = function(...){
      grDevices::dev.off()
      if (private$whether_contour) {
        private$plot_predict_contour(data = self$testing, ...)
      }else{
        private$plot_predict_line(...)
      }
      invisible(self)
    },
    plot_silent = function(logic){
        private$silent <- logic
        invisible(self)
    },
    # Methods for SGSIM.
    normlize = function(){
      if (!private$whether_predict) {
        cat("GPR: Please predict first")
      }
      sd <- self$para["sd_t"]
      self$normalized_n <- private$make_norm(dat = self$n_sws, 
                                             trend = self$testing, sd0 = sd)
      private$plot_norm()
      invisible(self)
    },
    make_sgsim_file = function(flname, dat = self$normalized_n, dig = 6,
                               fltital = "kaminokoike"){
      infoout <- data.frame(X=c(fltital, 3, "x"), Z=c("","","z"), N=c("","","N"))
      dataout <- data.frame(X = dat$x, z = dat$z, N = round(dat$nsws, digits = dig))
      write.table(infoout,file = flname,sep = "           ",
                  col.names=FALSE,row.names = FALSE,quote = FALSE)
      write.table(dataout,file = flname,sep = "           ",
                  col.names=FALSE,row.names = FALSE,quote = FALSE,append = TRUE)
    }
  ),
  #
  private = list(
    kernel_fun = "kernel_g",
    raw_n_sws = NULL,
    limit_x = c(-1,200),
    mesh_size_v = 0,
    mesh_size_h = 0,
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
    whether_opt = 0, 
    whether_predict = FALSE,
    whether_rigid = FALSE,
    silent = FALSE,
    ######################### Pure Functions ##############################
    # Read data from a file. 
    read_MAIC = function(file_name){
      rawdat<-read.csv(file = file_name) 
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
    # Make rigid data.
    rigidify = function(dat = self$n_sws, size_v, from_zero){
      # Rigidify the vertical axis ("z") of input data.
      if(size_v <= 0){
        size_v <- dplyr::group_by(dat, x) |> 
          summarise(dif = z %>% diff() %>% abs() %>% mean()) |> 
          dplyr::select(dif) |> unlist() |> mean()
        dplyr::ungroup(dat)
      }
      dat <- dplyr::group_by(dat, x) |> 
        dplyr::mutate(z = size_v*(1:length(x))) 
      if (from_zero) {
        dat <- dat |> dplyr::mutate(z = z - size_v)
      }
      dplyr::ungroup(dat)
      return(dat)
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
      if (size <= 0) {
        points_new <- points
      }else{
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
      }
      return(points_new)
    },
    #
    lerp = function(x, y, new_x){
      # Linear interpolation
      interp_f <- approxfun(x, y)
      new_y <- interp_f(new_x)
      return(new_y)
    },
    find_n = function(df, x, z){
      # Find N value from x and z coordinates.
      indx <- which(df$x == x & df$z == z)
      return(df$nsws[indx])
    },
    ####################### End of pure functions #########################
    # 
    mesh_data = function(size_v, size_h){
   # Prepare testing data (mesh).
      depth_raw <- self$n_sws |> 
        dplyr::group_by(x) |> 
        dplyr::summarise(min_depth = min(z), max_depth = max(z)) |>
        dplyr::ungroup() |>
        dplyr::arrange(x)
      # Calculate the unobserved points according to the "size_h" argument.
      x_denser <- private$divide_KP(points = depth_raw$x, size = size_h)
      depth <- data.frame(x = x_denser,
                          min = private$lerp(depth_raw$x, depth_raw$min_depth, new_x = x_denser),
                          max = private$lerp(depth_raw$x, depth_raw$max_depth, new_x = x_denser))
      # Determine the default type of result figure.
      if (nrow(depth) > nrow(depth_raw)) {
        private$whether_contour <- TRUE
      }else{
        private$whether_contour <- FALSE
      }
      #  Calculate the mean of interval on the vertical. 
      #  It will be used for unobserved points when users do not assign a vertical mesh size.
      if (size_v <= 0 & size_h > 0) {
        size_v_ex <- self$n_sws |> dplyr::group_by(x) |> 
          dplyr::summarise(dif= diff(z) %>% mean()) |>
          dplyr::ungroup() |>
          dplyr::select(dif) |>
          unlist() |> mean()
      }else{
        size_v_ex <- size_v
      }
      # Calculate mesh of unobserved points
      mesh_ex <- dplyr::filter(depth,!(x %in% depth_raw$x)) |> 
        as.list() |>
        rlang::set_names(NULL) |>
        purrr::pmap_dfr(function(dis, min, max){
          value <- seq(min, max, by = size_v_ex)
          len <- length(value)
          tibble::tibble(x = rep(dis,times = len), z = value)
        }) |>
        dplyr::mutate(y=0)
      # Create mesh of observed points.
      mesh_raw <- list(x= depth_raw$x, kp_z = split(self$n_sws$z, self$n_sws$x)) |> 
        rlang::set_names(NULL) |>
        purrr::pmap_dfr(function(dis, kp){
          kp <- unlist(kp)
          value <- private$divide_KP(points = kp, size = size_v)
          len <- length(value)
          tibble::tibble(x = rep(dis, times = len), z = value)
        })|>
        dplyr::mutate(y=0)
      # bind the 2 kinds of meshes.
      testing_mesh <- dplyr::bind_rows(mesh_raw, mesh_ex) |> 
        dplyr::arrange(x,z)
      #
      return(testing_mesh)
    },
    #
    check_predict = function(){
     if(purrr::some(self$para, is.na)){
       stop("Error in GPR: pleas use 'set_parameter' method to set ALL parameters (set nu only when you want to use WM model).")
     }
     if(!private$whether_set_parameter){
       stop("Error in GPR: it seems that the opt method has been interrupted for some reason. Please set the parameters or optimize them once more.")
     }
     if(private$kernel_fun == "kernel_wm"){
       if (private$whether_nu  == FALSE) {
         stop("Error in GPR: you have not set the initial value of parameter nu.")
       }
     }else{
       if (private$whether_nu  == TRUE) {
         warning("Warning in GPR: you set a value of nu without using WM model.")
       }
     }
     if(!private$whether_mesh){
       stop("Error in GPR: please use 'create_mesh' to creat a mesh before predicting.")
     }
    },
    # ------ plotting functions ------
    #
    plot_raw = function(){
      self$raw_pic = private$plot_points(dat = self$n_sws)
      if (!private$silent) {
        print(self$raw_pic)
      }
    },
    plot_norm = function(){
      self$norm_pic = private$plot_points(dat = self$normalized_n)
      if (!private$silent) {
        print(self$norm_pic)
      }
    },
    plot_points = function(dat){
      fig = ggplot2::ggplot(data = dat)+
        ggplot2::geom_point(mapping = aes(x=nsws,y=z))+
        ggplot2::geom_line(mapping = aes(x=nsws,y=z),orientation = "y")+
        ggplot2::scale_y_reverse()+
        ggplot2::facet_wrap(~x, nrow = 2)
    },
    #
    plot_predict_line = function(nrow = 2){
      self$test_pic <- ggplot() +
        geom_point(data = self$n_sws, mapping = aes(x = nsws, y = z))+
        geom_line(data = self$testing, mapping = aes(x = nsws, y = z),color="#D3323F",orientation = "y")+
        scale_y_reverse()+
        facet_wrap(~x, nrow = nrow)
      #
      if (!private$silent) {
        self$test_pic |> print()
      }
    },
    #
    plot_predict_contour = function(data, TitleColorbar, TitleX, TitleY, begin_X, max_Z,
                                    font_family = 'Times New Roman',font_size = 30, width = 1200 , height = 600){
      #
      DAT <- interp::interp(x = data$x, y = data$z, z = data$nsws, 
                            nx = max(self$n_sws$x), ny = 40, method="linear")
      DAT$z <- t(DAT$z)
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
      Y_min<-min(DAT$y)
      DAT$y<- DAT$y-Y_min
      #DAT$y<-rev(DAT$y)   #reverse the value of Y axes
      Y_max<-max(DAT$y)
      # shift the X axes
      shif_X <- fig_arguments$shifX
      DAT$x<-DAT$x+shif_X
      X_min<-min(DAT$x)
      # range of Z axes
      if(missing(max_Z)) max_Z <- max(DAT$z)
      Z_max <- max_Z
      #
      fig <- plot_ly(
        type = 'contour',
        z = DAT$z,
        y = DAT$y,
        x = DAT$x,
        autocontour = F,
        colorscale= list(c(0, 0.33,0.66, 1), c('#FF0000', '#FFFB00', '#339502', '#023AB2')),
        contours = list(
          start = 0,
          end = Z_max,
          size = Z_max/10,
          showlines = FALSE
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
        dtick = 10,                 
        ticklen = 8,              
        tickfont = font_for_tick,   
        tickwidth = 1,              
        tickcolor = toRGB("black"),
        zeroline = FALSE
      )
      y <- list(
        title=fig_arguments$title_Y,
        titlefont = list(size=font_size),  
        ticklen = 5,
        tickfont = font_for_tick,  
        tickwidth = 1,
        range = c(Y_max,0)
      )
      self$test_pic <- fig %>% layout(xaxis = x, yaxis = y,
                                      font = list(family = font_family ,size= font_size))
      #
      if (!private$silent) {
        self$test_pic |> print()
      }
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
      k11 <- `diag<-`(k11, diag(k11) + k11[1,1]*0.2)
      f <- -0.5 * t(private$z) %*% MASS::ginv(k11*k11_r_pivot) %*% private$z - 
        0.5 *(log(k11_r_pivot) * nrow(k11)+ log(det(k11))) + private$thirdterm
      # f <- -f
      f <- as.numeric(f)
      return(f)
    },
    #
    opt_bfgs = function(par, func, lower, upper){
      if (!private$whether_set_parameter) {
        stop("Error in GPR: pleas set initial parameters before using L-BFGS-B. ($set_parameter)")
      }
      if (!private$whether_set_scope) {
        stop("Error in GPR: pleas set the scope of each parameter before using L-BFGS-B. ($set_par_scope)")
      }
      func <- func |> private$more_info()
      try(self$result_BFGS <- optim(par, func, method = "L-BFGS-B", lower = lower, upper = upper))
    },
    #
    opt_ga = function(lower, upper, func){
      if (!private$whether_set_scope) {
        stop("Error in GPR: pleas set the scope of each parameter before using GA. ($set_par_scope)")
      }
      private$whether_set_parameter <- FALSE # If GA was interrupted, users can not predict with the initial parameters.
      self$result_ga <- GA::ga(type = "real-valued", fitness = func, 
                        lower = lower, upper = upper,
                        popSize = 120, maxiter = 100, run = 20, parallel = TRUE,
                        optim = FALSE)
      name_para <- names(self$para)
      self$para <- self$result_ga@solution[1,] |> as.numeric()
      names(self$para) <- name_para
      private$whether_set_parameter <- TRUE
    },
    #
    more_info = function(func){
      # A function operator used for likelihood function. 
      # Write the argument into 'self$para' and record the time that the function has been run.
      # The value will be reversed for optimizing.
      force(func)
      time <- 0
      function(x, ...){
        res <- func(x, ...)
        names(res) <- NULL
        self$para <- x
        time <<- time + 1
        cat("L-BFGS-B | iter = ", time, " | parameters: ", "\n" ,sep = "")
        print(x)
        cat("\n")
        -res      # BFGS method finds a minimum, but we want the maximum of likelihood.
      }
    },
    repdict_log = function(){
      cat("GPR: Prediction has been done. you can check the result through '$testing'.", "\n")
      cat("     Kernel function: ", private$kernel_fun, "\n", sep = "")
      cat("     Parameters: ", "\n", sep = "")
      print(self$para)
    },
   make_norm = function(dat, trend, sd0){
     len_dat <- unlist(dat$z) |> length()
     # compute the normalized N_sws
     for (i in seq_len(len_dat)){
       one_sample <- private$find_n(df = trend, x = dat$x[i], z = dat$z[i])
       dat$nsws[i] <- abs(one_sample - dat$nsws[i])/sd0
     }
     return(dat)
   }
  )
)
