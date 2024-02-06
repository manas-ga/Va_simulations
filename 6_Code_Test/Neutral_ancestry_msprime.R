###############################################################################################
######### Function to simulate neutral ancestry and simulate mutations using msprime ##########
###############################################################################################

# source("/mnt/c/Users/msamant/Documents/GitHub/Va_simulations/6_Code_test/Neutral_ancestry_msprime.R")
# python_path = "/home/manas_ga/miniconda3/envs/Va/bin/python"


### The function takes the following inputs:
# 1. sequence_length, 2. mu_msp (mutation rate for the msprime sim), 3. r_msp (recombination rate), 4. n_samples (number of diplod individuals to be sampled), 5. python_path (The path where python is tored on the system)


### Produces the following outputs
# 1. c (the genotype matrix with rows being haplotypes, and columns being loci), 2. pos (positions of mutations)

# For now I am assuming that the order of the loci in c is the same as the order in pos (to be confirmed !!!!!!!!)

###########################################################

# Install package "reticulate"

###########################################################

sim_msprime_ancecstry = function(sequence_length, 
                                 Ne, 
                                 mu_msp, 
                                 r_msp, 
                                 n_sample, # "n_samples" is the number of diploid individuals to be simulated (!Not genomes!)
                                 python_path) { # The path where python is stored
  
  
  # Specify the path to python
  reticulate::use_python(python_path)
  
  # Import msprime
  msprime <- reticulate::import("msprime")
  
  # Simulate the ancestry
  ts <- msprime$sim_ancestry(n_sample, sequence_length=sequence_length, recombination_rate=r_msp, population_size=Ne)
  
  # Simulate mutations
  mts = msprime$sim_mutations(ts, rate=mu_msp)
  
  print(paste("There are ", mts$num_sites, " sites with ", mts$num_mutations, " mutations"))
  
  # Store the c matrix
  c = t(mts$genotype_matrix()%%2)

  # Store the position
  pos = mts$sites_position
  
  # Output c matrix and positions
  return(list(c = c, pos = pos))
}

