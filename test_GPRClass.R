source("GPR_class.R")
test <- GPR$new("kaminokoike_SWS.dat", plot_silent = FALSE)$
  create_mesh()$
  set_par_scope(sof_h_t = c(0.1,100), sof_v_t = c(0.1, 10), sd_t = c(1,50),
                sof_h_r = c(0.1,100), sof_v_r = c(0.1, 2), sd_r = c(1,50))$
  opt("GA")$
  predict()

test$create_mesh(size_h = 10, size_v = 0)
test$set_parameter(sof_h_t = 32, sof_h_r = 0.5, sd_t = 5, sd_r = 5, sof_v_t = 2.6, sof_v_r = 0.2)
test$set_par_scope(sof_h_t = c(0.1,100), sof_v_t = c(0.1, 10), sd_t = c(1,50),
                    sof_h_r = c(0.1,100), sof_v_r = c(0.1, 2), sd_r = c(1,50))
test$opt("GA")
test$set_plot_lines(T)
test$predict()
test$plot_predict()
test$set_limit_x(c(-1,50))
save.image("1.RData")
