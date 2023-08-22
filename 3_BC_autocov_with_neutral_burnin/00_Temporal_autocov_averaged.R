# source("C:/Academics/Post-doc/Va_simulations/3_With_neutral_burnin_msprime/00_Temporal_autocov_averaged.R")


##########################################################################################################
####### Compute temporal autocovariances, averaged over many replicate simulations #######################
##########################################################################################################

## What does this script do?
# 1. Generate neutral diversity and attache selection coefficie s using "0_neutral_burnin.py"
# 2. Run a forward simulation using "1_FWD_simulation.slim"
# 3. Extract mutations using "2_Extract_mutations_genomes.py"
# 4. Calculate the temporal autocovariance btween start_gen and subsequent generations using the function  "compute_temp_autocov()"
# stored in "C:/Academics/Post-doc/Va_simulations/3_With_neutral_burnin_msprime/a_Functions/Function_to_compute_temporal_autocovariances.R"
# 5. Repeat 1 to 4 many times and average the temporal autocovariances calculated in each replicate




########################################################################
########### paths of various scripts and functions #####################
########################################################################

msprime_path = "C:/Academics/Post-doc/Va_simulations/3_With_neutral_burnin_msprime/0_neutral_burnin.py"

slim_path = "C:/Academics/Post-doc/Va_simulations/3_With_neutral_burnin_msprime/1_FWD_simulation.slim"

extract_mut_path = "C:/Academics/Post-doc/Va_simulations/3_With_neutral_burnin_msprime/2_Extract_mutations_genomes.py"

source("C:/Academics/Post-doc/Va_simulations/3_With_neutral_burnin_msprime/a_Functions/Function_to_compute_temporal_autocovariances.R")

mutations_path = "C:/Academics/Post-doc/Va_simulations/3_With_neutral_burnin_msprime/b_Interim_files/Mutations"

output_path = "C:/Academics/Post-doc/Va_simulations/3_With_neutral_burnin_msprime/c_Output"

##########################################################################################################################################
# Information on the SLiM simulation requited for the "compute_temp_autocov()" function

pop_size = 1000
start_gen = 1
end_gen = 51



# Create an empty vector temporal autocovariances
# In each replicate the temporal autocovariance vector for that replicate will be added using cbind()
# After all the replicates of done average temporal autocovariances will be calculated using rowMeans on this matrix

tac_matrix = c()




###################################################################################
############################ Loop over replicates #################################
###################################################################################


# Number of replicates

n_reps = 30

for (rep in 1:n_reps){
  
  message(paste("Replicate ", rep, " in progress..."))
    
  #######################
  ##### run msprime #####
  #######################  
    
  ## msprime assigns seed randomly, so no need to worry
  
  message("Genrating neutral diversity using msprime...")
    
  system(paste("python", msprime_path))
  
  #######################
  ###### run SLiM #######
  ####################### 
  
  # SLiM by default generates random seeds ina very odd way (based on the clock and process ID)
  # Generate a random seede for SLiM
  
  slim_seed = sample(1:100000, 1)
  
  message("Running the forward simulation using SLiM...")
  
  system(paste(paste("slim -seed", slim_seed), slim_path))
  
  
  ###############################
  ###### extract mutations ######
  ###############################
  
  message("extracting mutations...")
  
  system(paste("python", extract_mut_path))
  
  
  ##########################################
  ### Calculate temporal autocovariances ###
  ##########################################
  
  message("calculating temporal autocovariance...")
  
  tac = compute_temp_autocov(start_gen, end_gen, mutations_path, pop_size)
  tac_matrix = cbind(tac_matrix, tac)
  
  }

# Calculate average temporal autocovariance across replicates ignoring NAs

tac_av = rowMeans(tac_matrix, na.rm=TRUE)



pdf(paste(output_path, "/Average_temp_cov_s_0.05_g_1e-09.pdf", sep = ""), onefile = F)
plot(tac_av, xlab = "Generation", ylab = paste("Temporal autocovariance averaged over ", n_reps, " replicates"))
abline(0,0)
dev.off()

plot(tac_av, xlab = "Generation", ylab = paste("Temporal autocovariance averaged over ", n_reps, " replicates"))
abline(0,0)