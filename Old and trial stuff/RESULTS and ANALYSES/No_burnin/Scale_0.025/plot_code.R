library(ggplot2)

#### Code for making plots

d = read.csv("Data_no_burnin.csv", header=T)

d$n_cages = factor(d$n_cages, levels = c("3", "5", "10"))
d$n_ind_exp = factor(d$n_ind_exp, levels = c("100", "500", "1000"))
d$ngen_expt = factor(d$ngen_expt, levels = c("1", "3", "5"))



d_f = d[d$bdelta_method=="fixed",]
d_e = d[d$bdelta_method=="estimate",]


##################################################################################################
####### Standard plots (r_expt = 1.4e-6, n_ind_expt = 1000, n_cages = 10, ngen_expt = 3) #########
##################################################################################################

# bdelta fixed at c(0,0)

d_f_std = d_f[d_f$ngen_expt=="3"&d_f$n_ind_exp=="1000"&d_f$n_cages=="10"&d_f$r_expt==1.4e-06,]
p_f_std = ggplot(d_f_std, aes(y = vA_est, x = vA_true))
p_f_std + theme_bw() + geom_point() + geom_abline(intercept = 0, slope = 1) + labs(x = "True Vw", y = "Estimate of Vw") + theme(text = element_text(size = 15))

# bdelta to be estimated

d_e_std = d_e[d_e$ngen_expt=="3"&d_e$n_ind_exp=="1000"&d_e$n_cages=="10"&d_e$r_expt==1.4e-06,]
p_e_std = ggplot(d_e_std, aes(y = vA_est, x = vA_true))
p_e_std + theme_bw() + geom_point() + geom_abline(intercept = 0, slope = 1) + labs(x = "True Vw", y = "Estimate of Vw") + theme(text = element_text(size = 15))

##################################################################################################


#################################################
########### A. bdelta fixed at c(0,0) ###########
#################################################

### Map length ###

d_f_r = d_f[d_f$ngen_expt=="3"&d_f$n_ind_exp=="1000"&d_f$n_cages=="10",]
p_f_r = ggplot(d_f_r, aes(y = vA_est, x = vA_true, color = as.character(r_expt*sequence_length)))
p_f_r + theme_bw() + geom_point() + geom_abline(intercept = 0, slope = 1) + labs(x = "True Vw", y = "Estimate of Vw", color = "Map length (M)") + scale_color_manual(values = c("#0072B2", "#999999", "#009E73")) + theme(text = element_text(size = 15))

### n_ind_expt ###

d_f_ind = d_f[d_f$ngen_expt=="3"&d_f$n_cages=="10"&d_f$r_expt==1.4e-06,]
p_f_ind = ggplot(d_f_ind, aes(y = vA_est, x = vA_true, color = n_ind_exp))
p_f_ind + theme_bw() + geom_point() + geom_abline(intercept = 0, slope = 1) + labs(x = "True Vw", y = "Estimate of Vw", color = "Population size") + scale_color_manual(values = c("#0072B2", "#999999", "#009E73")) + theme(text = element_text(size = 15))

### n_cages ###

d_f_cage = d_f[d_f$ngen_expt=="3"&d_f$n_ind_exp=="1000"&d_f$r_expt==1.4e-06,]
p_f_cage = ggplot(d_f_cage, aes(y = vA_est, x = vA_true, color = n_cages))
p_f_cage + theme_bw() + geom_point() + geom_abline(intercept = 0, slope = 1) + labs(x = "True Vw", y = "Estimate of Vw", color = "Replicate populations") + scale_color_manual(values = c("#0072B2", "#999999", "#009E73")) + theme(text = element_text(size = 15))

### ngen_expt ###

d_f_gen = d_f[d_f$n_ind_exp=="1000"&d_f$n_cages=="10"&d_f$r_expt==1.4e-06,]
p_f_gen = ggplot(d_f_gen, aes(y = vA_est, x = vA_true, color = ngen_expt))
p_f_gen + theme_bw() + geom_point() + geom_abline(intercept = 0, slope = 1) + labs(x = "True Vw", y = "Estimate of Vw", color = "Generations") + scale_color_manual(values = c("#0072B2", "#999999", "#009E73")) + theme(text = element_text(size = 15))


#################################################
########### B. bdelta to be estimated ###########
#################################################

### Map length ###

d_e_r = d_e[d_e$ngen_expt=="3"&d_e$n_ind_exp=="1000"&d_e$n_cages=="10",]
p_e_r = ggplot(d_e_r, aes(y = vA_est, x = vA_true, color = as.character(r_expt*sequence_length)))
p_e_r + theme_bw() + geom_point() + geom_abline(intercept = 0, slope = 1) + labs(x = "True Vw", y = "Estimate of Vw", color = "Map length (M)") + scale_color_manual(values = c("#0072B2", "#999999", "#009E73")) + theme(text = element_text(size = 15))

### n_ind_expt ###

d_e_ind = d_e[d_e$ngen_expt=="3"&d_e$n_cages=="10"&d_e$r_expt==1.4e-06,]
p_e_ind = ggplot(d_e_ind, aes(y = vA_est, x = vA_true, color = n_ind_exp))
p_e_ind + theme_bw() + geom_point() + geom_abline(intercept = 0, slope = 1) + labs(x = "True Vw", y = "Estimate of Vw", color = "Population size") + scale_color_manual(values = c("#0072B2", "#999999", "#009E73")) + theme(text = element_text(size = 15))

### n_cages ###

d_e_cage = d_e[d_e$ngen_expt=="3"&d_e$n_ind_exp=="1000"&d_e$r_expt==1.4e-06,]
p_e_cage = ggplot(d_e_cage, aes(y = vA_est, x = vA_true, color = n_cages))
p_e_cage + theme_bw() + geom_point() + geom_abline(intercept = 0, slope = 1) + labs(x = "True Vw", y = "Estimate of Vw", color = "Replicate populations") + scale_color_manual(values = c("#0072B2", "#999999", "#009E73")) + theme(text = element_text(size = 15))

### ngen_expt ###

d_e_gen = d_e[d_e$n_ind_exp=="1000"&d_e$n_cages=="10"&d_e$r_expt==1.4e-06,]
p_e_gen = ggplot(d_e_gen, aes(y = vA_est, x = vA_true, color = ngen_expt))
p_e_gen + theme_bw() + geom_point() + geom_abline(intercept = 0, slope = 1) + labs(x = "True Vw", y = "Estimate of Vw", color = "Generations") + scale_color_manual(values = c("#0072B2", "#999999", "#009E73")) + theme(text = element_text(size = 15))

