
# Rscript /mnt/c/Users/msamant/Documents/GitHub/Va_simulations/1_Simulations/Misc_scripts_files/Create_parameter_grid_std.R

########################################################################
##### Script to write a text file containing the grid of parameters ####
########################################################################

# Col 1 = mu
# Col 2 = r*sequence_length
# Col 3 = r_expt*sequence_length
# Col 4 = n_ind_exp
# Col 5 = n_cages
# Col 6 = ngen2
# Col 7 = flip_sel_coef
# Col 8 = mut_ratio

nsims = 50 # number of simulations for each set
mu_list = seq(3e-8, 2e-7, length = nsims)

param_matrix = data.frame(matrix(NA, nrow = nsims, ncol  = 8))

param_matrix[,1] = mu_list
param_matrix[,2] = 1.4
param_matrix[,3] = 1.4
param_matrix[,4] = 1000
param_matrix[,5] = 10
param_matrix[,6] = 3
param_matrix[,7] = 0
param_matrix[,8] = 0


write.table(noquote(param_matrix), file = "/mnt/c/Users/msamant/Documents/GitHub/Va_simulations/1_Simulations/000_parameter_grid.txt", sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)
