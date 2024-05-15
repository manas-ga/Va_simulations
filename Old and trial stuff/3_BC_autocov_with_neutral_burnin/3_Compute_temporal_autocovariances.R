#########################################################
################### Specify paths #######################
#########################################################

# Working directory

setwd("C:/Academics/Post-doc/Va_simulations/3_With_neutral_burnin_msprime")

# Add the source of the file containing the function to compute va, q1, q2, q3, etc.

source("a_Functions/Function_to_Compute_V1_Q1_Q2_Q3.R")

# Specify the folder containing c matrices for genomes

genomes_path = "b_Interim_files/C_Matrices"

# Specify the folder containing mutations files

mutations_path = "b_Interim_files/Mutations"

####################################################################################
##### Input information regarding generations for which data has been recorded #####
####################################################################################


start_gen = 1
end_gen = 51

# Population Size

pop_size = 1000




#################################################
####### Compute temporal autocovariances ########
#################################################

### Create an empty matrix to store allelic frequencies at neutral loci in each generation
# Only include neutral loci segregating in generation 1
# If some of the loci get fixed/lost in subsequent generations, NAs should be inserted

P = c()

# Read the mutations in start_gen

mut_0 = read.csv(paste(mutations_path, "/mutations_", start_gen, ".txt", sep = ""), sep = " ")

# Filter out selected loci

mut_0_neutral = mut_0[mut_0$s==0,]

# store the frequencies of neutral mutations in start_gen in P
# Frequency = (Number of genomes)/(2*popsize)

P = cbind(P, mut_0_neutral$Number/(2*pop_size))



### Loop through the remaining files identifying neutral mutations from start_gen and recording their frequencies

for (gen in (start_gen+1):end_gen){
  # Read the file storing mutation information
  mut = read.csv(paste(mutations_path ,"/mutations_", gen, ".txt", sep = ""), sep = " ")
  
  # Create an empty vector to store frequencies of mutations in the current generation
  freq = c()
  
  # Loop through the permanent IDs of neutral mutations segregating in start_gen
  # i.e. Loop through Permanent IDs in mut_0_neutral
  # Check if each mutation is present in the current generation
  # If present record the frequency in freq, otherwise add NA to freq
  
  for(mutation in mut_0_neutral$Permanent_ID){
    if(mutation %in% mut$Permanent_ID){freq = rbind(freq, (mut[which(mut$Permanent_ID==mutation),]$Number)/(2*pop_size))}
    else {freq = rbind(freq, NA)}
  }
  
  # Add the vector freq to P
  
  P = cbind(P, freq)
  print(paste(round((gen-start_gen)*100/(end_gen-start_gen), 2), " % complete!"))
  
}


# Calculate the matrix of allele frequency changes

del_P = P[,2:(end_gen - start_gen + 1)] - P[,1:(end_gen - start_gen)] 


# Calculate the covariance matrix for allele frequency change

Q = cov(del_P, use = "pairwise.complete")


# Plot the temporal autocovariance of generation 1 with all other generations

#pdf("c_Output/Q_0.pdf", onefile = F)
plot(Q[1,]/(mean(P[,1])*(1-mean(P[,1]))),
     xlab = "Generation", ylab = "Temporal autocovariance with the first generation") 
#dev.off()


# Plot the temporal autocovariance of the middle generation with subsequent generations
# Q is scaled by the diversity at the initial generation

#pdf("c_Output/Q_mid.pdf", onefile = F)
plot(Q[round(end_gen/2, 0),]/(mean(complete.cases(P[,round(end_gen/2, 0)]))*(1-mean(complete.cases(P[,round(end_gen/2, 0)])))), 
     xlim = c(round(end_gen/2, 0), end_gen), 
     xlab = "Generation", ylab = paste("Temporal autocovariance with generation ", round(end_gen/2, 0))) 
#dev.off()



### BC (2020) say in the Methods that temporal autocovariances can be generated simply as a byproduct of tracking the minor allele frequency
### Therefore, it is more appropriate to randomly assign the reference allele
### Therefore, multiplying randomly chosen half the rows of del_P by -1 and recalculating temporal autocovariances

# Randomly select rows to be multiplied by -1

#rand_rows = sample(1:nrow(P), nrow(P)/2, replace = F)

#del_P_corr = del_P
#del_P_corr[rand_rows,] = del_P_corr[rand_rows,]*(-1)

#Q_corr = cov(del_P_corr, use = "pairwise.complete")

#plot(Q_corr[1,]/(mean(P[,1])*(1-mean(P[,1])))) # Q_corr scaled by the diversity in start_gen


#plot((Q_corr[1,]/(mean(P[,1])*(1-mean(P[,1])))) , (Q[1,]/(mean(P[,1])*(1-mean(P[,1])))))


#### Plot the number of segregating loci for each generation by counting NAs in the columns of P
# Create an empty vector to store this info

seg_sites = c()

# Loop through columns of P, counting NAs and building up seg_sites

for (gen in 1:ncol(P)){seg_sites[gen] = nrow(P) - sum(is.na(P[, gen]))}


#pdf("c_Output/segregating_sites.pdf", onefile = F)
plot(seg_sites, ylab = "Number of neutral segregating sites", xlab = "Generation")
#dev.off()



#######################################################################################
############### Compute initial Va (for a sample of 1000 individuals) #################
#######################################################################################

# Calculates 1. Va, 2. Q1, 3. Q2, 4. Q3, 5. summed site diversity, 6. no. of segregating loci, 
# 7. no. of selected loci, 7. variance in alphas, 8. mean of alphas

#BC_quantities = compute_va_q1_q2_q3(paste(genomes_path, "/c_matrix_", start_gen, ".csv", sep = ""), 
#                                    paste(mutations_path, "/mutations_", start_gen, ".txt", sep = ""))


# Calculate additive genic variance in start_gen

mut_0$Frequency = mut_0$Number/(2*pop_size)
mut_0$Diversity = 2*(mut_0$Frequency)*(1-mut_0$Frequency)
va_0 = sum((mut_0$Diversity)*((mut_0$s)/2)^2)
    
