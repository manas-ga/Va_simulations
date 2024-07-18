########################################################################
##### Script to write a text file containing the grid of parameters ####
########################################################################

# Col 1 = mu
# Col 2 = r*sequence_length
# Col 3 = r_expt*sequence_length
# Col 4 = n_ind_exp
# Col 5 = n_cages
# Col 6 = ngen_expt

nsims = 10 # number of simulations for each set
mu_list = seq(3.45e-07, 3.45e-06, length = nsims)

param_matrix = matrix(NA, nrow = nsims, ncol  = 7)

param_matrix[,1] = mu_list
param_matrix[,2] = 1.4
param_matrix[,3] = 1.4
param_matrix[,4] = 1000
param_matrix[,5] = 10
param_matrix[,6] = 3
param_matrix[,7] = "fixed"


write.table(param_matrix, file = "/mnt/c/Users/msamant/Documents/GitHub/Va_simulations/5_History_sim/000_parameter_grid.txt", sep = " ", col.names = FALSE, row.names = FALSE)
