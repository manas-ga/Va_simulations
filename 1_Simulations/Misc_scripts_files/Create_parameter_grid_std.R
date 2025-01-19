
# Rscript /mnt/c/Users/msamant/Documents/GitHub/Va_simulations/1_Simulations/Misc_scripts_files/Create_parameter_grid_std.R

########################################################################
##### Script to write a text file containing the grid of parameters ####
########################################################################

# Col 1 = mu_msp
# Col 2 = r*sequence_length
# Col 3 = r_expt*sequence_length
# Col 4 = n_ind_exp
# Col 5 = n_cages
# Col 6 = ngen2
# Col 7 = flip_sel_coef
# Col 8 = mut_ratio

standard_only = TRUE

## The standard set (to be used as a reference to compare)

nsims = 8 # number of simulations for each set
end_gen = 25000

mu_msp_list = if(end_gen==2){seq(3e-9, 2.35e-8, length = nsims)}else{seq(3.6e-8, 3.6e-7, length = nsims)}
mut_ratio = if(end_gen==2){1}else{0.0000}


param_matrix = data.frame(matrix(NA, nrow = nsims, ncol  = 8))

param_matrix[,1] = mu_msp_list
param_matrix[,2] = 250
param_matrix[,3] = 2
param_matrix[,4] = 1000
param_matrix[,5] = 10
param_matrix[,6] = 4
param_matrix[,7] = 0
param_matrix[,8] = mut_ratio

# Add to this matrix if simulations other than the standard set are also required

if(!standard_only){
  
  # Vary map_length (sequence_length*r)
  
  param_matrix_ml_v1 = param_matrix
  param_matrix_ml_v1[,2] = 50
  
  param_matrix_ml_v2 = param_matrix
  param_matrix_ml_v2[,2] = 100

  
  if(TRUE==FALSE){
    
    # Vary map_length_expt (sequence_length*r_expt)
    
    param_matrix_ml_v1 = param_matrix
    param_matrix_ml_v1[,3] = 0.01
    
    param_matrix_ml_v2 = param_matrix
    param_matrix_ml_v2[,3] = 0.2
    
    # Vary n_ind_exp
    
    param_matrix_nind_100 = param_matrix
    param_matrix_nind_100[,4] = 100
    
    param_matrix_nind_500 = param_matrix
    param_matrix_nind_500[,4] = 500
    
    # Vary n_cages
    
    param_matrix_cage_3 = param_matrix
    param_matrix_cage_3[,5] = 3
    
    param_matrix_cage_5 = param_matrix
    param_matrix_cage_5[,5] = 5
    
    # Vary ngen2
    
    param_matrix_ngen2_2 = param_matrix
    param_matrix_ngen2_2[,6] = 2
    
    param_matrix_ngen2_6 = param_matrix
    param_matrix_ngen2_6[,6] = 6
  
  
  
  # Combine all these matrices one below the other
  
  param_matrix = rbind(param_matrix,
                       param_matrix_ml_v1, param_matrix_ml_v2,
                       param_matrix_nind_100, param_matrix_nind_500,
                       param_matrix_cage_3, param_matrix_cage_5,
                       param_matrix_ngen2_2, param_matrix_ngen2_6)
  
  }
  param_matrix = rbind(param_matrix, param_matrix_ml_v1, param_matrix_ml_v2)
  
}

write.table(noquote(param_matrix), file = "/mnt/c/Users/msamant/Documents/GitHub/Va_simulations/1_Simulations/000_parameter_grid.txt", sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)
