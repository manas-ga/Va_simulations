
extract_slim_data = function(Set_ID,                # The unique ID of the set of simulations that are controlled by a single R script
                             sim = 1,               # Each set can have multiple sims, but - on the cluster sim must always 1
                             ngen2_optional = NULL, # Allows del_P to be calculated between ngen1 and manually specified ngen2 (which can be different from the last generation)
                             unzip = FALSE,         # Should the SLiM output file be unzipped, read, and then zipped back?
                             slim_output_path,      # The directory where the SLiM outputs (for parents and experimental replicates) are stored (as .txt files)
                             sim_param_path,        # The path to the directory where the .csv file containing simulation parameters is stored
                             extract_genomes_path,  # The path to the python script that extracts genomes and mutations from SLim outputs (3_Extract_genomes.py)
                             extract_mut_path,      # The path to the python script that extracts mutations from SLim outputs (2_Extract_mulations.py)
                             mutations_path,        # The directory where extracted mutations are to be stored (temp files)
                             c_matrix_path,         # The directory where extracted genomes are to be stored (temp files)
                             n_sample=NULL,         # Number of individuals sampled from the parents' generation (useful if n_ind_exp is large)
                             randomise = TRUE,      # Optionally the reference allele can be randomised
                             delete_temp_files = TRUE, 
                             pool_seq = FALSE, # Should the function simulate_pool_seq be used to sample allele frequencies in the experiment?
                             read_length = NULL,
                             coverage = NULL,
                             V_logmean = NULL,
                             verbose=TRUE
                             ){     
  
  ### Checks
  
  if(pool_seq){
    if(is.null(read_length)|is.null(coverage)|is.null(V_logmean)){stop("If poos_seq is true, read_length, coverage, and V_logmean must be provided")}
  }
  
  ################################
  ## Read simulation parameters ##
  ################################
  
  sim_params = read.csv(paste(sim_param_path, "/", Set_ID, "_Data.csv", sep = ""), header = T)
  # Restrict the data to current sim if the data frame has data from multiple sims
  if("sim"%in%colnames(sim_params)){sim_params = sim_params[sim_params$sim==sim,]}
  
  if(nrow(sim_params)==0){stop("No data found!")}
  
  # Load sim parameters
  
  n_ind_exp = sim_params$n_ind_exp     # Number of individuals in the experiment
  n_cages = sim_params$n_cages         # Number of replicates in the experiment
  if("ngen1"%in%colnames(sim_params)){ngen1 = sim_params$ngen1}else{ngen1 = 1}
  if("ngen2"%in%colnames(sim_params)){ngen2 = sim_params$ngen2}else{ngen2 = sim_params$ngen_expt + 1}
  ngen_expt = ngen2 - ngen1     # Number of generations over which allele frequency changes are calculated
  end_gen = sim_params$end_gen         # The generation number of the parents' generation
  
  if(!is.null(ngen2_optional)){
    if(ngen2_optional<=ngen1 | ngen2_optional>ngen2){stop("ngen2_optional must be greater than ngen1 and less than or equal to ngen2")}
    message(paste("Manually setting ngen2 to", ngen2_optional, "..."))
    ngen2 = ngen2_optional
    sim_params$ngen2 = ngen2 # Alsi change ngen2 in the sim_params
  }
  
  
  # If n_sample is not provided extract the genomes of all individuals in the parents' generation
  
  if(is.null(n_sample)){n_sample = n_ind_exp}
  if(n_sample>n_ind_exp){stop("n_sample cannot be greater than n_ind_exp")}
  
  if(verbose){
    message("Reading the state of the population in the parent's generation...")
    message("Extracting mutations and genomes from the parents' generation...")
  }
  
  # Unzip the SLiM output file
  
  if(unzip){
    if(verbose)message("Unzipping the SLiM output file...")
    system(paste("gunzip", paste(slim_output_path, "/", Set_ID, "_sim", sim, "_output_parents.txt.gz", sep = "")))
  }
  
  system(paste("python", 
               extract_genomes_path,                        # Path of the python script (3_Extract_genomes.py)
               paste(slim_output_path, "/", Set_ID, "_sim", sim, "_output_parents.txt", sep = ""),  # Path of the .txt file containing the SLiM output for the parent's generation 
               paste(mutations_path, "/", Set_ID, "_sim", sim, "_mutations_parents.txt", sep = ""), # Path of the .txt output file containing the mutations in the parents' generation 
               paste(c_matrix_path,"/", Set_ID, "_sim", sim, "_c_matrix_parents.csv", sep = ""),    # Path of the .csv output file containing the c matrix for genomes in the parents' generation
               n_sample))                                                                 # Number of individuals to be sampled randomly to construct the c matrix (just for space issues in case n_ind_exp is very large). Typically should be set to same as n_ind_exp
  
  # Rezip the SLiM output file
  
  if(unzip){
    if(verbose)message("Re-zipping the unzipped SLiM output file...")
    system(paste("gzip", paste(slim_output_path, "/", Set_ID, "_sim", sim, "_output_parents.txt", sep = "")))
  }
  
  
  # Read genomes
  c_genome = read.csv(paste(c_matrix_path,  "/", Set_ID, "_sim", sim, "_c_matrix_parents.csv", sep =""), header=F) # as.integer done to avoid scientific notation
  
  # Delete the .csv file
  if(delete_temp_files){system(paste("rm", paste(c_matrix_path,  "/", Set_ID, "_sim", sim, "_c_matrix_parents.csv", sep =""), sep = " "))}
  
  c_genome = as.matrix(c_genome)
  
  # Convert genome data into individual data (rows are individuals and columns are allele counts {0,1 or 2} at various sites). 
  # Note that genome 1 and genome 2 are from individual 1; 3 and 4 are from individual 2, and so on
  
  n0_individuals = nrow(c_genome)/2
  
  if(n0_individuals!=n_ind_exp){stop("The number of individuals obtained from c_genome do not match with n_ind_exp as logged during the simulation")}
  
  # If one samples individuals from the parents' generation while building the c matrix (i.e. when sample_size is less than n_ind_exp), the sample may not contain some low frequency mutations, i.e. some loci are not segregating in the sample, but are in the parents' population
  
  # identify the loci that are missing in the sample
  
  retained_loci = which(colSums(c_genome)!=0)
  
  # Trim the c matrix to contain only the retained loci, i.e. the loci that are segregating in the sample
  
  c_genome = c_genome[,retained_loci]
  n_sites = ncol(c_genome)
  
  ## Read the list of mutations output generated by SLiM and subsequently cropped out as a separate file by Python (for the parents' generation)
  
  pbar0 = colMeans(c_genome)
  # vector of starting allele frequencies
  
  mutations_0 = read.table(paste(mutations_path, "/", Set_ID, "_sim", sim, "_mutations_parents.txt", sep = ""), header=T, sep = " ")
  
  # Delete the .txt file
  if(delete_temp_files){system(paste("rm ", mutations_path, "/", Set_ID, "_sim", sim, "_mutations_parents.txt", sep = ""))}
  
  # Sort mutations based on the Temp_ID, so that the order of loci matches the order in c_ind
  
  mutations_0 = mutations_0[order(mutations_0$Temp_ID),]
  
  # Trim mutations to contain only the retained loci, i.e. the loci that are segregating in the sample
  
  mutations_0 = mutations_0[retained_loci,]
  
  # Check if allele frequencies calculated from c_genome match allele frequencies in mutations_0
  
  if(sum(pbar0 != mutations_0$Number/(2*n0_individuals)) == 0){message("Allele frequencies in c_genome match with allele frequencies from mutations_0")}else{stop("Allele frequencies in c_genome do not match with allele frequencies from mutations_0")}
  
  
  list_alpha = 2*(mutations_0$s)   # Vector of alphas
  SNPs = mutations_0$Position      # Vector of positions of mutations in the parents
  
  ############################################################################
  ############################################################################
  ####### Calculate allele frequencies in the experiment for each cage #######
  ############################################################################
  ############################################################################  
  
  # Create empty vector to create the data frame containing the following variables as columns:
  # 1. Raw delta P
  # 2. Projected delta P
  # 3. Cage ID (replicate)
  # 4. Locus ID
  
  P_matrix = c()  
  
  for (cage in (1:n_cages)){
    
    ###############################
    ###### extract mutations ######
    ###############################
    
    
    # There are two command line arguments (1. Path of the SLiM output file, 2. Path where mutations are to be written) 
    
    if(verbose){
      message(paste("Extracting mutations and storing allele frequencies in cage ", cage, " of simulation ", sim, "...", sep = ""))
    }
    
    for (gen in (end_gen + 1):(end_gen + ngen2)){
      
      # Unzip SLiM output files
      
      if(unzip){
        if(verbose){message("Unzipping the SLiM output file...")}
        system(paste("gunzip", paste(slim_output_path, "/", Set_ID, "_sim", sim, "_cage", cage, "_output_experiment_", as.integer(gen), ".txt", sep = "")))
      }
      
      # If pool_seq = TRUE extract both mutations and genomes, otherwise only extract mutations
      if(pool_seq){
        system(paste("python", 
                     extract_genomes_path,                        # Path of the python script (3_Extract_genomes.py)
                     paste(slim_output_path, "/", Set_ID, "_sim", sim, "_cage", cage, "_output_experiment_", as.integer(gen), ".txt", sep = ""),  # Path of the .txt file containing the SLiM output for the current cage and the current generation 
                     paste(mutations_path, "/", Set_ID, "_sim", sim, "_cage", cage, "_mutations_", as.integer(gen), ".txt", sep = ""), # Path of the .txt output file containing the mutations for the current cage and the current generation 
                     paste(c_matrix_path,"/", Set_ID, "_sim", sim, "_cage", cage, "_c_matrix_", as.integer(gen), ".csv", sep = ""),    # Path of the .csv output file containing the c matrix for genomes for the current cage and the current generation
                     n_sample)) 
        
      }else{
        system(paste("python", 
                     extract_mut_path, 
                     paste(slim_output_path, "/", Set_ID, "_sim", sim, "_cage", cage, "_output_experiment_", as.integer(gen), ".txt", sep = ""), 
                     paste(mutations_path, "/", Set_ID, "_sim", sim, "_cage", cage, "_mutations_", as.integer(gen), ".txt", sep = ""))) # as.integer done to avoid scientific notation
      }
      
      # Re-zip unzipped SLiM files
      if(unzip){
        if(verbose){message("Re-zipping the unzipped SLiM output file...")}
        system(paste("gzip", paste(slim_output_path, "/", Set_ID, "_sim", sim, "_cage", cage, "_output_experiment_", as.integer(gen), ".txt", sep = "")))
      }
    }
    
    
    ### Create an empty matrix to store allelic frequencies in each generation
    
    # If some of the loci get fixed/lost in subsequent generations, NAs should be inserted
    
    # Create an empty vector to store allele frequencies
    
    P = c()
    
    # store the frequencies in the parent's generation in P
    # Frequency = (Number of genomes)/(2*popsize)
    
    P = cbind(P, mutations_0$Number/(2*n_ind_exp))
    
    ### Loop through the generations of the experiment identifying mutations that were present in the parents' generation (using permanent IDs) and recording their frequencies
    
    for (gen in (end_gen+1):(end_gen+ngen2)){
      
      # Read the file storing mutation information
      mut = read.csv(paste(mutations_path, "/", Set_ID, "_sim", sim, "_cage", cage, "_mutations_", as.integer(gen), ".txt", sep = ""), sep = " ") # as.integer gets rid of scientific notation
      
      # Delete the .txt file
      if(delete_temp_files){system(paste("rm ", mutations_path, "/", Set_ID, "_sim", sim, "_cage", cage, "_mutations_", as.integer(gen), ".txt", sep = ""))}
      
      # Sort mutations based on the Temp_ID, so that the order of loci matches the order in the C and L matrices for this gen
      
      mut = mut[order(mut$Temp_ID),]
      
      #####################################################################
      ### Optionally resample allele frequencies by simulating pool-seq ###
      #####################################################################
      
      if(pool_seq){
        message(paste("Simulating pool-seq on cage ", cage, ", generation ",  gen, "...", sep = ""))
        # Load the genomes of the current cage
        
        c_genome_cage_gen = read.csv(paste(c_matrix_path,"/", Set_ID, "_sim", sim, "_cage", cage, "_c_matrix_", as.integer(gen), ".csv", sep = ""), header = F)
        pool_seq_data = simulate_pool_seq(c_genome = c_genome_cage_gen,
                                       SNPs = mut$Position,      
                                       sequence_length = sim_params$sequence_length,
                                       read_length = read_length,
                                       coverage = coverage,
                                       V_logmean = V_logmean)
        
        # Delete the .csv file
        if(delete_temp_files){system(paste("rm ", c_matrix_path,"/", Set_ID, "_sim", sim, "_cage", cage, "_c_matrix_", as.integer(gen), ".csv", sep = ""))}
        
        # Modify the Numbers column in mut based on the sampled frequencies
        mut$Number = pool_seq_data$p*(2*n_ind_exp)
      }
      
      # Create an empty vector to store frequencies of mutations in the current generation
      freq = c()
      
      # Loop through the permanent IDs of  mutations segregating in end_gen (parents' generation)
      # i.e. Loop through Permanent IDs in mutations_0
      # Check if each mutation is present in the current generation
      # If present, record the frequency in freq, otherwise add either 0 or 1 to freq using the round() function
      
      for(mutation in mutations_0$Permanent_ID){
        if(mutation %in% mut$Permanent_ID){freq = c(freq, (mut[which(mut$Permanent_ID==mutation),]$Number)/(2*n_ind_exp))
        }else{
          freq = c(freq, round(P[which(mutations_0$Permanent_ID==mutation), gen - end_gen]))
          
        }
      }
      
      # Add the vector freq to P
      
      P = cbind(P, freq)
      
    }
    
    P_matrix = rbind(P_matrix, P)
    
  }
  
  
  ##### Calculate matrices of allele frequencies with rows as replicates ######
  
  pbar1 = matrix(NA, nrow = n_cages, ncol = n_sites)
  pbar2 = matrix(NA, nrow = n_cages, ncol = n_sites)
  for (i in 1:n_cages){
    pbar1[i,] = P_matrix[((i-1)*n_sites + 1):((i-1)*n_sites + n_sites), (ngen1 + 1)] ## Matrix of frequencies at the start of the experiment (ie F1 generation)
  }
  
  for (i in 1:n_cages){
    pbar2[i,] = P_matrix[((i-1)*n_sites + 1):((i-1)*n_sites + n_sites), (ngen2 + 1)] ## Matrix of frequencies at the end of the experiment
    
    
  }
  
  if(randomise){
    
    #########################################################
    ######## Randomise the reference allele in c_ind ########
    #########################################################
    
    # Randomly change the reference allele
    # This can be done as follows:
    # Add 0s to the allele counts of those alleles that stay the same, and -2 to those alleles that are to be switched
    # Then take a mod
    # Remember c_ind is a matrix of 0s 1s and 2s
    
    if(verbose){
      message("Randomising reference alleles in c_ind...")
    }
    
    # Generate a random vector of 0s (for no change) and -1s (for loci where the reference allele is to be switched)
    
    ran_vect = sample(c(0, -1), ncol(c_genome),  replace = T) 
    
    # Create a matrix with with as many rows as c_ind. 
    # Each row of this matrix should be made up of two times ran_vect (since we are working with allele counts, not frequencies). 
    # Because the same changes need to be applied to each individual
    
    ran_matrix = t(matrix(ran_vect, nrow = ncol(c_genome), ncol = nrow(c_genome)))
    
    # Calculate the allele counts of the new (randomised) reference alleles
    
    c_genome = abs(c_genome + ran_matrix)
    
    ##############################################################
    ######## Randomise the reference allele in list_alpha ########
    ##############################################################
    
    if(verbose){
      message("Randomising reference alleles in list_alpha...")
    }
    
    # Create a new vector of -1s (wherever ran_vect has -1) and 1s (wherever ran_vect has 0)
    ran_vect_alpha = ifelse(ran_vect==-1, -1, 1)
    
    list_alpha = list_alpha*ran_vect_alpha
    
    
    ################################################################
    ######## Randomise the reference allele in pbar1 and pbar2 #####
    ################################################################
    
    # Randomly change the reference allele
    # This can be done as follows:
    # Add 0s to the frequencies of those alleles that stay the same, and -1 to those alleles that are to be switched
    # Then take a mod
    
    if(verbose){
      message("Randomising reference alleles in pbar1 and pbar2...")
    }
    
    # Create a matrix with with as many rows as pbar1. Each row of this matrix should be made up of ran_vect (computed while randomising c_ind). 
    # Because the same changes need to be applied to each cage
    
    ran_matrix_pbar = t(matrix(ran_vect, nrow = ncol(pbar1), ncol = nrow(pbar1)))
    
    # Calculate the frequencies of the new (randomised) reference alleles
    
    pbar0 = abs(pbar0 + ran_vect)     
    pbar1 = abs(pbar1 + ran_matrix_pbar)
    pbar2 = abs(pbar2 + ran_matrix_pbar)      
    
  }
  
  return(list(c_genome=c_genome, list_alpha=list_alpha, SNPs=SNPs, ngen1=ngen1, ngen2=ngen2,  pbar0=pbar0, pbar1=pbar1, pbar2=pbar2, sim_params=sim_params))
  
}