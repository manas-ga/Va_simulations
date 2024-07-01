########################################################################
##### Script to write a text file containing the grid of parameters ####
########################################################################

nsims = 100 # number of simulations for each set
mu_list = seq(2.25e-09, 2.25e-08, length = nsims)

param_matrix = matrix(NA, nrow = 9, ncol  = 4)

param_matrix[,1] = rep(1.4, 9)
param_matrix[2,1] = 0.14
param_matrix[3,1] = 0.014

param_matrix[,2] = rep(1000, 9)
param_matrix[4,2] = 500
param_matrix[5,2] = 100

param_matrix[,3] = rep(10, 9)
param_matrix[6,3] = 5
param_matrix[7,3] = 3

param_matrix[,4] = rep(3, 9)
param_matrix[8,4] = 1
param_matrix[9,4] = 5

# Repeat the matrix nsims times

param_matrix_full = c()

for (sim in 1:nsims){
  param_matrix_full = rbind(param_matrix_full,param_matrix)
}

mu_list = rep(mu_list, nrow(param_matrix))
mu_list = sort(mu_list)

param_matrix_full = cbind(mu_list, param_matrix_full)
param_matrix_full = param_matrix_full[order(param_matrix_full[,1]),]

param_matrix_full = param_matrix_full[seq(1, 9*nsims, 9),] ## TEMP trial selection (TO BE DELETED !!!!!)

#param_matrix = rbind(param_matrix, param_matrix)
#param_matrix = cbind(c(rep("fixed", 9), rep("estimate", 9)), param_matrix)



write.table(param_matrix_full, file = "/mnt/c/Users/msamant/Documents/GitHub/Va_simulations/5_History_sim/000_parameter_grid.txt", sep = " ", col.names = FALSE, row.names = FALSE)
