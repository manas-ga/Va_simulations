library(ggplot2)

#### Code for making plots

d = read.csv("Data_no_burnin.csv", header=T)

# Delete row names in every alternate line

d = d[grep("Set_ID", d$Set_ID, invert=TRUE),]

d$vA_est = as.numeric(d$vA_est)
d$vA_true = as.numeric(d$vA_true)
d$vA_left = as.numeric(d$vA_true)
d$pdelta_est = as.numeric(d$pdelta_est)
d$bdelta_intercept_est = as.numeric(d$bdelta_intercept_est)
d$bdelta_slope_est = as.numeric(d$bdelta_slope_est)
d$r_expt = as.numeric(d$r_expt)
d$sequence_length = as.numeric(d$sequence_length)
d$mu = as.numeric(d$mu)

d$n_cages = factor(d$n_cages, levels = c("3", "5", "10"))
d$n_ind_exp = factor(d$n_ind_exp, levels = c("100", "500", "1000"))
d$ngen_expt = factor(d$ngen_expt, levels = c("1", "3", "5"))





##################################################################################################
####### Standard plots (r_expt = 1.4e-6, n_indxpt = 1000, n_cages = 10, ngen_expt = 3) #########
##################################################################################################


d_std = d[d$ngen_expt=="3"&d$n_ind_exp=="1000"&d$n_cages=="10"&d$r_expt==1.4e-06,]
p_std = ggplot(d_std, aes(y = vA_est, x = vA_true))
p_std + theme_bw() + geom_point() + geom_abline(intercept = 0, slope = 1) + labs(x = "True Vw", y = "Estimate of Vw") + theme(text = element_text(size = 15))

##################################################################################################

### Map length ###

d_r = d[d$ngen_expt=="3"&d$n_ind_exp=="1000"&d$n_cages=="10",]
p_r = ggplot(d_r, aes(y = vA_est, x = vA_true, color = as.character(r_expt*sequence_length)))
p_r + theme_bw() + geom_point() + geom_abline(intercept = 0, slope = 1) + labs(x = "True Vw", y = "Estimate of Vw", color = "Map length (M)") + scale_color_manual(values = c("#0072B2", "#999999", "#009E73")) + theme(text = element_text(size = 15))

### n_ind_expt ###

d_ind = d[d$ngen_expt=="3"&d$n_cages=="10"&d$r_expt==1.4e-06,]
p_ind = ggplot(d_ind, aes(y = vA_est, x = vA_true, color = n_ind_exp))
p_ind + theme_bw() + geom_point() + geom_abline(intercept = 0, slope = 1) + labs(x = "True Vw", y = "Estimate of Vw", color = "Population size") + scale_color_manual(values = c("#0072B2", "#999999", "#009E73")) + theme(text = element_text(size = 15))

### n_cages ###

d_cage = d[d$ngen_expt=="3"&d$n_ind_exp=="1000"&d$r_expt==1.4e-06,]
p_cage = ggplot(d_cage, aes(y = vA_est, x = vA_true, color = n_cages))
p_cage + theme_bw() + geom_point() + geom_abline(intercept = 0, slope = 1) + labs(x = "True Vw", y = "Estimate of Vw", color = "Replicate populations") + scale_color_manual(values = c("#0072B2", "#999999", "#009E73")) + theme(text = element_text(size = 15))

### ngen_expt ###

d_gen = d[d$n_ind_exp=="1000"&d$n_cages=="10"&d$r_expt==1.4e-06,]
p_gen = ggplot(d_gen, aes(y = vA_est, x = vA_true, color = ngen_expt))
p_gen + theme_bw() + geom_point() + geom_abline(intercept = 0, slope = 1) + labs(x = "True Vw", y = "Estimate of Vw", color = "Generations") + scale_color_manual(values = c("#0072B2", "#999999", "#009E73")) + theme(text = element_text(size = 15))

