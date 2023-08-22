rm(list = ls())
# source("C:/Academics/Post-doc/Va_simulations/4_Va_theory/00_Va_theory.R")

####################################
######### Packages #################
####################################

library(matlib)
library(MCMCglmm)
library(asreml)
library(Matrix)

########################################################################
########### paths of various scripts and functions #####################
########################################################################

msprime_path = "C:/Academics/Post-doc/Va_simulations/4_Va_theory/0_neutral_burnin.py"

slim_path = "C:/Academics/Post-doc/Va_simulations/4_Va_theory/1_FWD_simulation.slim"

extract_mut_path = "C:/Academics/Post-doc/Va_simulations/4_Va_theory/2_Extract_mutations_genomes.py"

#source("C:/Academics/Post-doc/Va_simulations/4_Va_theory/a_Functions/Function_to_compute_temporal_autocovariances.R")

mutations_path = "C:/Academics/Post-doc/Va_simulations/4_Va_theory/b_Interim_files/Mutations"
genomes_path = "C:/Academics/Post-doc/Va_simulations/4_Va_theory/b_Interim_files/C_Matrices"

output_path = "C:/Academics/Post-doc/Va_simulations/4_Va_theory/c_Output"

##########################################################################################################################################
# Information on the SLiM simulation requited for the "compute_temp_autocov()" function

pop_size = 1000
start_gen = 1 # The parental generation (no selection)
end_gen = 5
r = 1e-09 # Recombination rate (per site per generation)
n_cages = 10
n_sims = 15
sinc<-1e-6            # increment added to diagonal of allelic effect covariance matrix to make pd
projection = "LoM" ## What projection to use (L or LoM)

## Create an empty vector to store va_estimated and va_true in every simulation

va_est = rep(NA, n_sims)
va_true = rep(NA, n_sims)

# Loop over replicate simulations, each time incrementing the selection coefficient assigned in msprime

for (sim in 1:n_sims){
  
  message("##################################################")
  message(paste("Simulation number ", sim, " in progress..."))
  message("##################################################")

  ##########################################################################################################################################
  
  #######################
  ##### run msprime #####
  #######################  
    
  ## msprime assigns seed randomly, so no need to worry
  
  message("Generating neutral diversity using msprime...")
  
  # Assign the selection coefficient (the first argument of the commant that calls the msprime python script)
  
  s = 0.3*(sqrt(sim)/n_sims)
  #s = 0.25
    
  system(paste("python", msprime_path, s))
  
  
  
  
  #######################
  ###### run SLiM #######
  ####################### 
  
  # SLiM by default generates random seeds i na very odd way (based on the clock and process ID)
  # Generate a random seed for SLiM
  
  slim_seed = sample(1:100000, 1)
  
  message("Running the forward simulation using SLiM...")
  
  system(paste(paste("slim -seed", slim_seed), slim_path))
  
  
  ###############################
  ###### extract mutations ######
  ###############################
  
  message("extracting mutations...")
  
  system(paste("python", extract_mut_path))
  
  ##########################################################################################################################################
  
  
  #########################################################################
  ##### Calculate C, L, non-recombinant fraction and Va for start_gen #####
  #########################################################################
  
  # Read genomes
  c_genome = read.csv(paste(genomes_path, "/c_matrix_1.csv", sep =""), header=F)
  
  c_genome = as.matrix(c_genome)
  
  #Convert genome data into individual data (rows are individuals and columns are allele counts {0,1 or 2} at various sites). 
  # Note that genome 1 and genome 2 are from individual 1, and so on
  
  message("Calculating the C matrix for the starting population...")
  
  n_individuals = nrow(c_genome)/2
  n_sites = ncol(c_genome)
  
  c_ind = c_genome[seq(1, 2*n_individuals, 2),] + c_genome[seq(2, 2*n_individuals, 2),]
  
  
  
  #### Calculate the matrix of second mixed moments of c, ie. the L matrix
  
  message("Calculating the L matrix for the starting population...")
  
  L = cov(c_ind/2)
  
  # Check for 0s on the diagonal of L
  #which(diag(L)==0)
  
  ####### Notice that L has some 0s on the diagonal
  # This is happening because when I sample individuals in SLiM from a population of 10000, the sample does not contain some low frequency mutations, ie. some loci are not segregating in the sample, but are in the original population
  
  # identify the loci that are missing in the sample
  
  missing_loci = which(diag(L)==0)
  retained_loci = which(diag(L)!=0)
  
  # Trim L to contain only the retained loci
  
  L_ret = L[retained_loci,retained_loci]
  
  # Trim the c matrix to contain only the retained loci, i.e. the loci that are segregating in the sample
  
  c_ind_ret = c_ind[,retained_loci]
  
  # Correlation matrix
  
  #R = cov2cor(L_ret)
  
  
  ## Read the list of mutations output generated by SLiM and subsequently cropped out as a separate file by Python
  # From the first generation
  
  mutations_1 = read.table(paste(mutations_path, "/mutations_1.txt", sep = ""), header=T, sep = " ")
  
  # Sort mutations based on the Temp_ID, so that the order of loci matches the order in the C and L matrices
  
  mutations_1 = mutations_1[order(mutations_1$Temp_ID),]
  
  # Trim mutations to contain only the retained loci, i.e. the loci that are segregating in the sample
  
  mutations_ret = mutations_1[retained_loci,]
  
  # Calculate the number of loci that are retained
  
  n_sites_ret = ncol(c_ind_ret)
  
  ######## Eigen decomposition of the projection matrix
  # Can project on L or LoM depending on the value of the parameter projection
  
  
  if (projection == "L"){ ## Project on L
    message("Performing eigen-decomposition of L...")
    eigen_L = eigen(L)
  }
  
  else{ ## Project on LoM
  
    ### Construct the recombination map matrix (NRF)
    # Diagonal elements are 1
    # Off-diagonal elements are non-recombinant fractions
    
    # Construct the  non-recombination fraction matrix
    # Calculating distances between sites, multiply by the recombination rate
    # Subtract from 1
    
    message("Calculating the the matrix of non-recombinant fractions for the starting population...")
    
    NRF = c()
    
    # Loop through columns 
    
    for (site in 1:nrow(L)){
      dist = abs(mutations_1$Position - mutations_1$Position[site])
      rf = 0.5*(1 - exp(-2*dist*r))
      NRF = cbind(NRF, 1 - rf)
    }
    
    # Calculate the covariance structure for del_P under drift and recombibation 
    # var(del_P) = L/2N*SUM(1 to n_gen){(NRF^(n+1)(1 - 1/2N)^n)}
    
    LoM = matrix(0, nrow(L), nrow(L))
    
    for (n in 1:(end_gen - (start_gen + 1))){
      LoM = LoM + (1/(2*pop_size))*L*(NRF^(1+n))*(1 - 1/(2*pop_size))^n
      #LoR = LoR + L*(NRF^(1+n))*(1 - 1/(2*pop_size))^n
    }
    
    # Eigen-decomposition of L
    
    message("Performing eigen-decomposition of LoM...")
    
    eigen_L = eigen(LoM)
  
  }
  
  ### Eigen decomposition is expressed in our theory as LoR = U%*%D^2%*%t(U)
  # U is the matrix of eigenvecors
  # D is a diagonal matrix of the square roots of the eigenvalues
  # Retain only those eigen values (and corresponding eigen vectors) greater than tol
  
  tol  = sqrt(.Machine$double.eps)
  
  eigen_ret = which(eigen_L$values>tol) 
  
  U = eigen_L$vectors
  U = U[, eigen_ret]
  D = diag(sqrt(eigen_L$values[eigen_ret]))
  inv_D = diag(1/(diag(D)))
  
  
  ### Calculate fitnesses of individuals 
  
  message("Computing the actual Va for the starting population...")
  
  # Create an empty vector of individual fitnesses with each element initialized to 1 
  #This will have to be 0 when using an additive model for fitness
  
  #w = rep(1, n_individuals)
  
  # Calculate the fitness of each individual (i.e. populate w) by multiplying the contributions of each locus
  # The product will have to be replaced by a sum when things are additive
  # k loops over individuals, and m loops over sites
  
  #for (k in 1:n_individuals){
    # loop through loci, summing fitness contributions
    #for (m in 1:n_sites_ret){
      # Calculate the fitness of each individual using a multiplicative model for fitness
      #w[k] = w[k] * (1 + ((2-c_ind_ret[k,m])*mutations_ret$h[m] + (c_ind_ret[k,m]-1)/2)*c_ind_ret[k,m]*mutations_ret$s[m]) 
      
    #}
    # Add a progress bar
    #print(paste("Calculating fitnesses: ", round(k*100/n_individuals, 2), "% complete"))
  #}
  
  ######### Alternative way of calulating fitnesses (Courtesy Jarrod!!)
  
  w<-apply(1 + (t(t(2-c_ind_ret)*mutations_ret$h)+(c_ind_ret-1)/2)*t(t(c_ind_ret)*mutations_ret$s), 1, prod)
  
  ### Perform linear regression to compute th average effects, alphas
  
  
  #Calculate relative fitness
  
  Fitness = w/mean(w)
  
  fit_1 <- lm(Fitness ~ ., data = data.frame(c_ind_ret/2)) # Dividing by 2 because the c matrix is in the form of 0s, 1s or 2s
  
  # Remove the first coefficient (i.e. the intercept), and store the alphas
  
  list_alpha = coef(fit_1)[-1] 
  
  # Since (typically) n_individuals < n_loci, many coefficients cannot be estimated, and list_alpha, etc have NAs
  
  # Convert NAs to 0s
  
  list_alpha[is.na(list_alpha)]<-0
  
  
  # Calculate Va
  
  # While calculating Va NAs in list_alpha were converted to 0s
  
  va = t(list_alpha)%*%L_ret%*%list_alpha
  
  va_true[sim] = va
  
  print("###############")
  print(va)
  print("###############")

  
  ################################################
  ##### Calculate P and del_P for each cage ######
  ################################################
  
  message("Computing allele frequency changes for each cage...")
  
  # Create empty vector to store projected allele frequency change
  
  del_P = c()
  del_P1 = c()
  
  for (cage in 1:n_cages){
    
    message(paste("Cage", cage, "of simulation", sim, "in progress..."))
    
    #######################
    ###### run SLiM #######
    ####################### 
    
    # SLiM by default generates random seeds in a very odd way (based on the clock and process ID)
    # Generate a random seede for SLiM
    
    slim_seed1 = sample(1:100000, 1)
    
    message("Running the forward simulation using SLiM...")
    
    system(paste(paste("slim -seed", slim_seed1), slim_path))
    
    
    ###############################
    ###### extract mutations ######
    ###############################
    
    message("extracting mutations...")
    
    system(paste("python", extract_mut_path))
    
    # Read the mutations file for the first generation (start_gen)
    
    mutations_1 = read.table(paste(mutations_path, "/mutations_1.txt", sep = ""), header=T, sep = " ")
    
    # Sort mutations based on the Temp_ID, so that the order of loci matches the order in the C and L matrices
    
    mutations_1 = mutations_1[order(mutations_1$Temp_ID),]
    
    ### Create an empty matrix to store allelic frequencies in each generation
    
    # If some of the loci get fixed/lost in subsequent generations, NAs should be inserted
    
    message("Calculating allele frequency change...")
    
    P = c()
    
    
    # store the frequencies in start_gen in P
    # Frequency = (Number of genomes)/(2*popsize)
    
    P = cbind(P, mutations_1$Number/(2*pop_size))
    
    
    
    ### Loop through the remaining files identifying mutations from start_gen and recording their frequencies
    
    for (gen in (start_gen+1):end_gen){
      # Read the file storing mutation information
      mut = read.csv(paste(mutations_path ,"/mutations_", gen, ".txt", sep = ""), sep = " ")
      
      # Sort mutations based on the Temp_ID, so that the order of loci matches the order in the C and L matrices
      
      mut = mut[order(mutations_1$Temp_ID),]
      
      # Create an empty vector to store frequencies of mutations in the current generation
      freq = c()
      
      # Loop through the permanent IDs of  mutations segregating in start_gen
      # i.e. Loop through Permanent IDs in mut_0
      # Check if each mutation is present in the current generation
      # If present, record the frequency in freq, otherwise add NA to freq
      
      for(mutation in mutations_1$Permanent_ID){
        if(mutation %in% mut$Permanent_ID){freq = rbind(freq, (mut[which(mut$Permanent_ID==mutation),]$Number)/(2*pop_size))}
        else {
          freq = rbind(freq, round(P[which(mutations_1$Permanent_ID==mutation), gen -1]))
          #print(round(P[which(mutations_1$Permanent_ID==mutation), gen -1]))
        }
      }
      
      # Add the vector freq to P
      
      P = cbind(P, freq)
      #print(paste(round((gen-start_gen)*100/(end_gen-start_gen), 2), " % complete!"))
      
    }
    #print("!!!!!!!!!!!!!!!!!!!! TEST !!!!!!!!!!!!!!!")
    #print(P[3,3])
    #print("!!!!!!!!!!!!!!!!!!!! TEST !!!!!!!!!!!!!!!")
    # Calculate the matrix of allele frequency changes between gen 2 and gen 3
    
    del_P_current = P[,end_gen] - P[,start_gen+1]
    
    del_P_current[is.na(del_P_current)]<-0
    
    
    
    
    # Calculate projected allele frequency change in the new basis defined by the eigen decomposition of LoR
    
    message("Calculating projected allele frequency change...")
    
    del_P_proj_current = inv_D%*%t(U)%*%del_P_current # inv_D = diag(1/(diag(D))). The fuction inv(D) would not work. Removing inv_D%*%
    
    del_P = rbind(del_P, cbind(del_P_proj_current, cage, 1:length(del_P_proj_current)))
    del_P1 = rbind(del_P1, cbind(del_P_current, cage, 1:length(del_P_current)))
    
  }
  
  # Create a data frame
  
  d = data.frame("del_P_proj" = del_P[,1], "Cage" = del_P[,2], "Locus" = del_P[,3])
  d1 = data.frame("del_P" = del_P1[,1], "Cage" = del_P1[,2], "Locus" = del_P1[,3])
  
  # Save the data frame
  
  write.csv(d, paste(output_path, "/del_P_data.csv", sep = ""))
  
  plot(d[d$Cage==1,]$del_P, d[d$Cage==2,]$del_P)
  abline(0,1)
  
  
  
  
  ##### Fit MCMCglmm model
  
  # Define the covariance structure for Locus
  
  if(projection == "L"){
        SC = D*D
  }
  
  else{
    SC = inv_D%*%t(U)%*%L%*%L%*%U%*%inv_D
  }
  
  invSC = solve(SC+diag(length(eigen_ret))*sinc) # Adding a small bit to the diagonal so that SC becomes invertible
  
  invSC<-as(invSC, "sparseMatrix")
  attr(invSC, "INVERSE")<-TRUE
  attr(invSC,"rowNames")<-as.character(eigen_ret) # Used for sparse-form matrices
  attr(invSC,"colNames")<-as.character(eigen_ret) # Used for sparse-form matrices
  dimnames(invSC) <- list(eigen_ret,eigen_ret)  # used for full-form matrices
  
  # Convert Locus to a factor
  
  d$Locus = factor(d$Locus)
  
  
  prior = list(R = list(V = 1, nu = 0), G = list(G1 = list(V = 1, nu = 1, alpha.mu=0, alpha.V=1000)))
  fit_mcmc = MCMCglmm(del_P_proj ~ 1, random = ~ Locus, ginverse=list(Locus=invSC), family="gaussian", prior = prior, data = d, pr = T) 
  
  va_est[sim] = mean(fit_mcmc$VCV[1])*sum(diag(L))/(end_gen - start_gen - 1)^2
  
  #m1<-asreml(del_P_proj~1, random = ~vm(Locus, invSC), data=d)
  #va_est[sim] = summary(m1)$varcomp[1,1]*sum(diag(L))/(end_gen - start_gen - 1)^2


}

plot(va_true, va_est)
abline(0, 1)


## Save

pdf(paste(output_path, "/Va_Output_", projection, ".pdf", sep = ""), onefile = F)
plot(va_true, va_est, xlab = "True value of Va", ylab = "Estimate of Va from the model")
abline(0, 1)
#dev.off()
#pred = predict(fit_mcmc, marginal = NULL)


# plot(d1[d1$Cage==7,]$del_P, d1[d1$Cage==3,]$del_P)