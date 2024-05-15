#########################################################
################### Specify paths #######################
#########################################################

# Add the source of the file containing the function to compute va, q1, q2, q3

source("a_Functions/Function_to_Compute_V1_Q1_Q2_Q3.R")

# Specify the folder containing c matrices for genomes

genomes_path = "b_Interim_files/C_Matrices"

# Specify the folder containing mutations files

mutations_path = "b_Interim_files/Mutations"

#######################################################################
### Build matching filepaths for mutations and c matrices           ###
### Calculate va, q1, q2, q3 for each gen and store in a data frame ###
#######################################################################

# Enter the starting generation, the ending generation and the frequency of the output in SLiM

start_gen = 101000
end_gen = 200000
output_freq = 1000

quantities_bc = data.frame()

for (gen in seq(start_gen, end_gen, output_freq)){
  
  c_matrix_filepath = paste(genomes_path, "/c_matrix_", format(gen, scientific = F), ".csv", sep = "")
  mutations_filepath = paste(mutations_path, "/mutations_", format(gen, scientific = F), ".txt", sep = "")
  # The gen inside the paste() function needs to be formatted because otherwise for some generations it prints the scientific format (eg 2e+05, etc.)
  
  quantities_bc = rbind(quantities_bc, c(gen, compute_va_q1_q2_q3(c_matrix_filepath, mutations_filepath)))

  # Just adding a progress bar
  progress = (gen - (start_gen-output_freq))/(end_gen - (start_gen-output_freq))*100
  print(paste(progress, "% complete!"))
  
  }

colnames(quantities_bc) = c("Generation", "Va", "Q1", "Q2", "Q3")

write.csv(quantities_bc, "c_Output/Quantities_BC.csv")
  

