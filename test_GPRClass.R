
test <- GPR$new("kaminokoike_SWS.dat")
test$set_parameter(sof_h_t = 30, sof_h_r = 0.5, sd_t = 2, sd_r = 5, sof_v_t = 5, sof_v_r = 0.1)
# test$set_par_scope(sof_h_t = c(1,30), sof_v_t = c(0.1, 10), sd_t = c(1,10),
#                    sof_h_r = c(1,30), sof_v_r = c(0.1, 10), sd_r = c(1,10))
test$create_mesh(size_h = 10)
test$predict()
test$opt(mode = "GA")
test$para

