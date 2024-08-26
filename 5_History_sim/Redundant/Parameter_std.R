
# Rscript /mnt/c/Users/msamant/Documents/GitHub/Va_simulations/5_History_sim/Redundant/Parameter_std.R

########################################################################
##### Script to write a text file containing the grid of parameters ####
########################################################################

# Col 1 = mu
# Col 2 = r*sequence_length
# Col 3 = r_expt*sequence_length
# Col 4 = n_ind_exp
# Col 5 = n_cages
# Col 6 = ngen_expt
# Col 7 = bdelta_method
# Col 8 = flip_sel_coef
# Col 9 = mut_ratio

nsims = 10 # number of simulations for each set
mu_list = seq(5.56e-07, 5.56e-06, length = nsims)

param_matrix = data.frame(matrix(NA, nrow = nsims, ncol  = 8))

param_matrix[,1] = mu_list
param_matrix[,2] = 1.4
param_matrix[,3] = 1.4
param_matrix[,4] = 1000
param_matrix[,5] = 10
param_matrix[,6] = 3
param_matrix[,7] = "estimate"
param_matrix[,8] = 0
param_matrix[,9] = 0.0002


write.table(noquote(param_matrix), file = "/mnt/c/Users/msamant/Documents/GitHub/Va_simulations/5_History_sim/000_parameter_grid.txt", sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)
