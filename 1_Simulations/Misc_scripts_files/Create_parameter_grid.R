
# Rscript /mnt/c/Users/msamant/Documents/GitHub/Va_simulations/1_Simulations/Misc_scripts_files/Create_parameter_grid.R

########################################################################
##### Script to write a text file containing the grid of parameters ####
########################################################################

# Col 1 = mu
# Col 2 = r*sequence_length
# Col 3 = r_expt*sequence_length
# Col 4 = n_ind_exp
# Col 5 = n_cages
# Col 6 = ngen_expt
# Col 7 = flip_sel_coef
# Col 8 = mut_ratio

nsims = 50 # number of simulations for each set
mu_list = seq(5.56e-07, 5.56e-06, length = nsims)
test = TRUE # If TRUE, only selects parameters for the "standard" simulation set

param_matrix = matrix(NA, nrow = 9, ncol  = 7)

# Col 2

param_matrix[,1] = rep(1.4, 9) 
param_matrix[2,1] = 0.14
param_matrix[3,1] = 0.014

# Col 3

param_matrix[,2] = param_matrix[,1] # For now recombination rates in both phases are the same!!!

# Col 4

param_matrix[,3] = rep(1000, 9)
param_matrix[4,3] = 500
param_matrix[5,3] = 100

# Col 5

param_matrix[,4] = rep(10, 9)
param_matrix[6,4] = 5
param_matrix[7,4] = 3

# Col 6

param_matrix[,5] = rep(3, 9)
param_matrix[8,5] = 1
param_matrix[9,5] = 5

# Col 7 

param_matrix[,6] = 0

# Col 8

param_matrix[,7] = 0

# Repeat the matrix nsims times

param_matrix_full = c()

for (sim in 1:nsims){
  param_matrix_full = rbind(param_matrix_full,param_matrix)
}

mu_list = rep(mu_list, nrow(param_matrix))
mu_list = sort(mu_list)

param_matrix_full = cbind(mu_list, param_matrix_full)
param_matrix_full = param_matrix_full[order(param_matrix_full[,1]),]

if(test){
  param_matrix_full = param_matrix_full[seq(1, 9*nsims, 9),] ## TEMP trial selection (TO BE DELETED !!!!!)
}

write.table(param_matrix_full, file = "/mnt/c/Users/msamant/Documents/GitHub/Va_simulations/1_Simulations/000_parameter_grid.txt", sep = " ", col.names = FALSE, row.names = FALSE)
