library(ggplot2)
library(cowplot)
library(latex2exp)

#### Code for making plots

d = read.csv("Data_r_12_July_2024.csv", header=T)

# Delete row names in every alternate line

d = d[grep("Set_ID", d$Set_ID, invert=TRUE),]

d$vA_est = as.numeric(d$vA_est)
d$vA_true = as.numeric(d$vA_true)
d$vA_left = as.numeric(d$vA_true)
d$pdelta_est = as.numeric(d$pdelta_est)
d$bdelta_intercept_est = as.numeric(d$bdelta_intercept_est)
d$bdelta_slope_est = as.numeric(d$bdelta_slope_est)
d$sigma2delta_est = as.numeric(d$sigma2delta_est)
d$r_expt = as.numeric(d$r_expt)
d$r = as.numeric(d$r)
d$sequence_length = as.numeric(d$sequence_length)
d$mu = as.numeric(d$mu)
d$s_pmq = as.numeric(d$s_pmq)
d$seg_sites = as.numeric(d$seg_sites)
d$seg_sites_neu = as.numeric(d$seg_sites_neu)
d$seg_sites_del = as.numeric(d$seg_sites_del)
d$seg_sites_ben = as.numeric(d$seg_sites_ben)
d$pdelta_emp = as.numeric(d$pdelta_emp)
d$bdelta_intercept_emp = as.numeric(d$bdelta_intercept_emp)
d$bdelta_slope_emp = as.numeric(d$bdelta_slope_emp)
d$sigma2delta_emp = as.numeric(d$sigma2delta_emp)
d$s_pmq = as.numeric(d$s_pmq)
d$n_cages = factor(d$n_cages, levels = c("3", "5", "10"))
d$n_ind_exp = factor(d$n_ind_exp, levels = c("100", "500", "1000"))
d$ngen_expt = factor(d$ngen_expt, levels = c("1", "3", "5"))

d$map_length_history = factor(d$r*d$sequence_length)
d$map_length_expt = factor(d$r_expt*d$sequence_length)

# V_A and seg sites

p_r_vA_true = ggplot(d, aes(y = vA_true, x = map_length_history)) + theme_bw() + geom_point() + labs(x = "Map length in the history phase (M)", y = TeX(r"(True $V_A$)")) + theme(text = element_text(size = 18))
p_r_vA_error = ggplot(d, aes(y = vA_est - vA_true, x = map_length_history)) + theme_bw() + geom_point() + labs(x = "Map length in the history phase (M)", y = TeX(r"(Error in the estimate of $V_A$)")) + theme(text = element_text(size = 18))
p_r_seg_sites = ggplot(d, aes(y = seg_sites, x = map_length_history)) + theme_bw() + geom_point() + labs(x = "Map length in the history phase (M)", y = "Number of segregating sites") + theme(text = element_text(size = 18))
p_r_seg_sites_del = ggplot(d, aes(y = seg_sites_del, x = map_length_history)) + theme_bw() + geom_point() + labs(x = "Map length in the history phase (M)", y = "Number of deleterious segregating sites") + theme(text = element_text(size = 18))

pdf("vA_seg_sites_combined.pdf", onefile = F, w=12, h=12)
plot_grid(p_r_seg_sites, p_r_seg_sites_del, p_r_vA_true, p_r_vA_error, labels = "AUTO")
dev.off()

# Empirical analysis parameters

p_r_pdelta_emp = ggplot(d, aes(y = pdelta_emp, x = map_length_history)) + theme_bw() + geom_point() + labs(x = "Map length in the history phase (M)", y = TeX(r"(Empirical $p_{\alpha} $)")) + theme(text = element_text(size = 18))

p_r_bdelta_intercept_emp = ggplot(d, aes(y = bdelta_intercept_emp, x = map_length_history)) + theme_bw() + geom_point() + labs(x = "Map length in the history phase (M)", y = TeX(r"(Empirical $\beta^{(0)}_{\bar{\alpha}}$)")) + theme(text = element_text(size = 18))

p_r_bdelta_slope_emp = ggplot(d, aes(y = bdelta_slope_emp, x = map_length_history)) + theme_bw() + geom_point() + labs(x = "Map length in the history phase (M)", y = TeX(r"(Empirical $\beta^{(1)}_{\bar{\alpha}}$)")) + theme(text = element_text(size = 18))

p_r_sigma2delta_emp = ggplot(d, aes(y = sigma2delta_emp, x = map_length_history)) + theme_bw() + geom_point() + labs(x = "Map length in the history phase (M)", y = TeX(r"(Empirical $\sigma^2_{\bar{\alpha}}$)")) + theme(text = element_text(size = 18))

pdf("Analysis_param_emp.pdf", onefile = F, h = 12, w = 12)
plot_grid(p_r_pdelta_emp, p_r_bdelta_intercept_emp, p_r_bdelta_slope_emp, p_r_sigma2delta_emp, labels = "AUTO")
dev.off()

# Estimate analysis parameters

p_r_pdelta_est = ggplot(d, aes(y = pdelta_est, x = map_length_history)) + theme_bw() + geom_point() + labs(x = "Map length in the history phase (M)", y = TeX(r"(Estimate of $p_{\alpha} $)")) + theme(text = element_text(size = 18))

p_r_bdelta_intercept_est = ggplot(d, aes(y = bdelta_intercept_est, x = map_length_history)) + theme_bw() + geom_point() + labs(x = "Map length in the history phase (M)", y = TeX(r"(Estimate of $\beta^{(0)}_{\bar{\alpha}}$)")) + theme(text = element_text(size = 18))

p_r_bdelta_slope_est = ggplot(d, aes(y = bdelta_slope_est, x = map_length_history)) + theme_bw() + geom_point() + labs(x = "Map length in the history phase (M)", y = TeX(r"(Estimate of $\beta^{(1)}_{\bar{\alpha}}$)")) + theme(text = element_text(size = 18))

p_r_sigma2delta_est = ggplot(d, aes(y = sigma2delta_est, x = map_length_history)) + theme_bw() + geom_point() + labs(x = "Map length in the history phase (M)", y = TeX(r"(Estimate of $\sigma^2_{\bar{\alpha}}$)")) + theme(text = element_text(size = 18))

pdf("Analysis_param_est.pdf", onefile = F, h = 12, w = 12)
plot_grid(p_r_pdelta_est, p_r_bdelta_intercept_est, p_r_bdelta_slope_est, p_r_sigma2delta_est, labels = "AUTO")
dev.off()

# Plotting empirical and estimated parameters side by side

d_r_pdelta = data.frame("map_length_history" = c(d$map_length_history, d$map_length_history), "pdelta" = c(d$pdelta_emp, d$pdelta_est), "Type" = c(rep("Empirical", nrow(d)), rep("Estimate", nrow(d))))
p_r_pdelta = ggplot(d_r_pdelta, aes(y=pdelta, x = map_length_history, color = Type)) + theme_bw() + geom_point(position=position_dodge(width=0.3)) + labs(x = "Map length in the history phase (M)", y = TeX(r"($p_{\alpha} $)"), color = "") + theme(text = element_text(size = 18)) 

d_r_bdelta_intercept = data.frame("map_length_history" = c(d$map_length_history, d$map_length_history), "bdelta_intercept" = c(d$bdelta_intercept_emp, d$bdelta_intercept_est), "Type" = c(rep("Empirical", nrow(d)), rep("Estimate", nrow(d))))
p_r_bdelta_intercept = ggplot(d_r_bdelta_intercept, aes(y=bdelta_intercept, x = map_length_history, color = Type)) + theme_bw() + geom_point(position=position_dodge(width=0.3)) + labs(x = "Map length in the history phase (M)", y = TeX(r"($\beta^{(0)}_{\bar{\alpha}}$)"), color = "") + theme(text = element_text(size = 18)) 

d_r_bdelta_slope = data.frame("map_length_history" = c(d$map_length_history, d$map_length_history), "bdelta_slope" = c(d$bdelta_slope_emp, d$bdelta_slope_est), "Type" = c(rep("Empirical", nrow(d)), rep("Estimate", nrow(d))))
p_r_bdelta_slope = ggplot(d_r_bdelta_slope, aes(y=bdelta_slope, x = map_length_history, color = Type)) + theme_bw() + geom_point(position=position_dodge(width=0.3)) + labs(x = "Map length in the history phase (M)", y = TeX(r"($\beta^{(1)}_{\bar{\alpha}}$)"), color = "") + theme(text = element_text(size = 18)) 

d_r_sigma2delta = data.frame("map_length_history" = c(d$map_length_history, d$map_length_history), "sigma2delta" = c(d$sigma2delta_emp, d$sigma2delta_est), "Type" = c(rep("Empirical", nrow(d)), rep("Estimate", nrow(d))))
p_r_sigma2delta = ggplot(d_r_sigma2delta, aes(y=sigma2delta, x = map_length_history, color = Type)) + theme_bw() + geom_point(position=position_dodge(width=0.3)) + labs(x = "Map length in the history phase (M)", y = TeX(r"($\sigma^2_{\bar{\alpha}}$)"), color = "") + theme(text = element_text(size = 18)) 


# Extract legend from one of the plots
legend = get_legend(p_r_pdelta + guides(color = guide_legend(nrow = 1)) + theme(legend.position = "bottom"))

# Make a combined plot without legends
plot_combined = plot_grid(p_r_pdelta + theme(legend.position="none"), 
                          p_r_bdelta_intercept + theme(legend.position="none"), 
                          p_r_bdelta_slope + theme(legend.position="none"), 
                          p_r_sigma2delta + theme(legend.position="none"), 
                          labels = "AUTO")

# Put together plot_combined and the legend using cowplot

pdf("Analysis_param_combined.pdf", onefile = F, h =12, w = 12)
plot_grid(plot_combined, legend, ncol = 1, rel_heights = c(1, .05))
dev.off()
