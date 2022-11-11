source("GPR_class.R")
神ノ子池 <- GPR$new("kaminokoike_SWS.dat")
神ノ子池$create_mesh()
神ノ子池$set_parameter(sof_h_t = 23.840455, sof_v_t = 6.672777, sd_t = 2, 
                sof_h_r = 72.307987, sof_v_r = 6.418074, sd_r = 2,
                nu = 5.949964)
神ノ子池$set_par_scope(sof_h_t = c(0.1,150), sof_v_t = c(0.1, 10), sd_t = c(1e-1,50),
                sof_h_r = c(0.1,100), sof_v_r = c(0.1, 10), sd_r = c(1e-1,50),
                nu = c(0.1,10))
神ノ子池$opt("BFGS")
神ノ子池$predict()

神ノ子池$create_mesh(size_h = 10, size_v = 0)
神ノ子池$set_parameter(sof_h_t = 23.840455, sof_v_t = 6.672777, sd_t = 4.8, 
                   sof_h_r = 72.307987, sof_v_r = 6.418074, sd_r = 4.801119,
                   nu = 5.949964)
# 神ノ子池$set_par_scope(sof_h_t = c(0.1,150), sof_v_t = c(0.1, 10), sd_t = c(1,50),
#                     sof_h_r = c(0.1,100), sof_v_r = c(0.1, 10), sd_r = c(1,50),
#                    nu = c(0.1,10))
# 神ノ子池$opt("GA")
神ノ子池$set_plot_lines(F)
神ノ子池$predict()
神ノ子池$plot_predict()
神ノ子池$set_limit_x(c(-1,50))
save.image("wm.RData")
