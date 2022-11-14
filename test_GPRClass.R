source("GPR_class.R")
神の子池 <- GPRSG$new("kaminokoike_SWS.dat")
神の子池$set_mesh_size(size_v = 0.25, size_h = 1)
神の子池$create_mesh(rigid = TRUE, from_zero = TRUE)
神の子池$set_parameter(sof_h_t = 23.840455, sof_v_t = 6.672777, sd_t = 2, 
                sof_h_r = 72.307987, sof_v_r = 6.418074, sd_r = 0.5,
                nu = 5.949964)
神の子池$set_par_scope(sof_h_t = c(0.1,150), sof_v_t = c(0.1, 10), sd_t = c(1e-1,50),
                sof_h_r = c(0.1,100), sof_v_r = c(0.1, 10), sd_r = c(1e-1,50),
                nu = c(0.1,10))
#神の子池$opt("GA")
神の子池$predict()
神の子池$normlize()
神の子池$make_sgsim_file(flname = "kaminoko_ike.dat", dig = 6)

save.image("wm.RData")

神の子池$set_sgs_mesh(nx = 117, ny = 35, dx = 1, dy = 0.25, num_sim = 60, 
                      members = 10)
神の子池$read_sgs(flname = "Sgsim result.dat")
神の子池$denormalize()
神の子池$take_exp()
神の子池$plot_sgs_sample(number =  1)
神の子池$plot_sgs_member(number =  10)




