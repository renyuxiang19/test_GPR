source("GPR_class.R")
test <- GPR$new("kaminokoike_SWS.dat")$
  create_mesh()$
  set_par_scope(sof_h_t = c(0.1,150), sof_v_t = c(0.1, 10), sd_t = c(1,50),
                sof_h_r = c(0.1,100), sof_v_r = c(0.1, 10), sd_r = c(1,50),
                nu = c(0.1,10))$
  opt("GA")$
  predict()

test$create_mesh(size_h = 10, size_v = 0)
test$set_parameter(sof_h_t = 23.840455, sof_v_t = 6.672777, sd_t = 4.149789, 
                   sof_h_r = 72.307987, sof_v_r = 6.418074, sd_r = 4.801119,
                   nu = 5.949964)
# test$set_par_scope(sof_h_t = c(0.1,150), sof_v_t = c(0.1, 10), sd_t = c(1,50),
#                     sof_h_r = c(0.1,100), sof_v_r = c(0.1, 10), sd_r = c(1,50),
#                    nu = c(0.1,10))
# test$opt("GA")
test$set_plot_lines(F)
test$predict()
test$plot_predict()
test$set_limit_x(c(-1,50))
save.image("wm.RData")
