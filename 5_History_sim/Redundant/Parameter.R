########################################################################
##### Script tp write a text file containing the grid of parameters ####
########################################################################


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

param_matrix[,4] = rep(9, 3)
param_matrix[8,4] = 1
param_matrix[9,4] = 5

param_matrix = rbind(param_matrix, param_matrix)

param_matrix = cbind(c(rep("fixed", 9), rep("estimate", 9)), param_matrix)



write.table(param_matrix, file = "param.txt", sep = " ", col.names = FALSE, row.names = FALSE)
