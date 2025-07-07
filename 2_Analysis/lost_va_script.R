# This script calculates:
# 1. True V_A and V_a using \alpha = \eta (old) or \alpha = \eta + 0.25(1 -2p)(\eta)^2 (new) 
# 2. Average va left (using old and new definitions of \alpha) at the end of the experiment
# 3. If finescale = TRUE, also records the \eta and frequency for each sefrefating locus.

rm(list = ls())
# Rscript /mnt/c/Users/msamant/Documents/GitHub/Va_simulations/2_Analysis/lost_va_script.R


##################
### File Paths ###
##################

if(Sys.info()["nodename"]!="SCE-BIO-C06645"){message("WARNING: Command line arguments required!!!")}

### On AC3 ###

# Old

if(Sys.info()["nodename"]%in%c("bigfoot", "bigshot", "bigbird", "bigyin", "biggar", "bigwig", "c1", "c2", "c3", "c4", "c5", "c6")){
  analysis_path = "/ceph/users/marun/Va_simulations/2_Analysis"
  file_storage_path = "/data/obbard/Va_simulations/analyses" # File storage path is the designated storage space on AC3 instead of the /home directory on qm
  Vw_library_path = "/ceph/users/marun/Va_simulations/Vw"
}

# New

if(Sys.info()["nodename"]%in%c("ac3-n1", "ac3-n2", "ac3-n3", "ac3-n4", "ac3-n5", "ac3-n6")){
  analysis_path = "~/Va_simulations/2_Analysis"
  file_storage_path = "/mnt/hel/obbard/Va_simulations/analyses" # File storage path is the designated storage space on AC3 instead of the /home directory on qm
  Vw_library_path = "~/Va_simulations/Vw"
}

### On Vera ###

if(Sys.info()["nodename"]=="vera.bio.ed.ac.uk"){
  analysis_path = "~/Va_simulations/2_Analysis"  
  file_storage_path = "/RawData/Manas_NERC_Simulations"
  Vw_library_path = "~/Va_simulations/Vw"
  
}

### On Eddie ###

if(grepl("ecdf.ed.ac.uk", Sys.info()["nodename"])){
  analysis_path = "/home/msamant/Va_simulations/2_Analysis"
  file_storage_path = "/exports/eddie/scratch/msamant"
  Vw_library_path = "/home/msamant/Va_simulations/Vw"
}

### On Manas's PC

if(Sys.info()["nodename"]=="SCE-BIO-C06645"){
  
  if(Sys.info()["sysname"]=="Linux"){   ## Local Wsl
    
    analysis_path = "/mnt/c/Users/msamant/Documents/GitHub/Va_simulations/2_Analysis" 
    file_storage_path = "/mnt/c/Users/msamant/Documents/GitHub/Va_simulations/1_Simulations"
    slim_path = "/mnt/u/Datastore/CSCE/biology/groups/hadfield/Va_simulations/sim_files"
    Vw_library_path = "/mnt/c/Users/msamant/Documents/GitHub/Va_simulations/Vw"
    
    
  }else{                                ## Local windows
    
    analysis_path = "C:/Users/msamant/Documents/GitHub/Va_simulations/2_Analysis" 
    file_storage_path = "C:/Users/msamant/Documents/GitHub/Va_simulations/1_Simulations"
    slim_path = "U:/Datastore/CSCE/biology/groups/hadfield/Va_simulations/sim_files"
    Vw_library_path = "C:/Users/msamant/Documents/GitHub/Va_simulations/Vw"
  }
  
}

### On Jarrod's PC ###

if(Sys.info()["nodename"]=="sce-bio-c04553"){  
  analysis_path = "~/Work/Va_simulations/2_Analysis"
  file_storage_path = "~/Work/Va_simulations/1_Simulations"
  slim_path = "/Volumes/hadfield/Va_simulations/sim_files"
  Vw_library_path = "~/Work/Va_simulations/Vw"
}

####################
### Script paths ###
####################

extract_genomes_path = file.path(analysis_path, "3_Extract_genomes.py")                                           ## Python script extracting mutations and genomes from the SLiM output file 
extract_mut_path = file.path(analysis_path, "2_Extract_mutations.py")                                             ## Python script extracting mutations from the SliM output file

temp_files_path = file.path(file_storage_path, "b_Interim_files")  
output_path = paste(file_storage_path, "/c_Output", sep = "")                                                     ## Path where .csv file(s) containing final data are to be stored 

# On Eddie use a different output path (since files in /scracth are deleted after a months)

if(grepl("ecdf.ed.ac.uk", Sys.info()["nodename"])){
  output_path = paste("~", "/c_Output", sep = "")  
}

###########################
### Path to SLiM output ###
###########################

test = TRUE

if(Sys.info()["nodename"]=="SCE-BIO-C06645"|Sys.info()["nodename"]=="sce-bio-c04553"){
  slim_output_path = if(test){file.path(file_storage_path, "b_Interim_files/SLiM_outputs")}else{file.path(slim_path, "SLiM_outputs")}
}

if(Sys.info()["nodename"]%in%c("bigfoot", "bigshot", "bigbird", "bigyin", "biggar", "bigwig", "c1", "c2", "c3", "c4", "c5", "c6", "ac3-n1", "ac3-n2", "ac3-n3", "ac3-n4", "ac3-n5", "ac3-n6")){
  slim_output_path = file.path(file_storage_path, "b_Interim_files/SLiM_outputs")
}

##########################################
### Path to SLiM simulation parameters ###
##########################################

if(Sys.info()["nodename"]=="SCE-BIO-C06645"|Sys.info()["nodename"]=="sce-bio-c04553"){
  sim_param_path = if(test){file.path(file_storage_path, "c_Output")}else{file.path(slim_path, "sim_params")}   
}

if(Sys.info()["nodename"]%in%c("bigfoot", "bigshot", "bigbird", "bigyin", "biggar", "bigwig", "c1", "c2", "c3", "c4", "c5", "c6", "ac3-n1", "ac3-n2", "ac3-n3", "ac3-n4", "ac3-n5", "ac3-n6")){
  sim_param_path = file.path(file_storage_path, "c_Output")
}

# Load packages and functions

library(RhpcBLASctl)

# Control the number of BLAS threads if running on a cluster
if(Sys.info()["nodename"]!="SCE-BIO-C06645"|Sys.info()["nodename"]!="sce-bio-c04553"){blas_set_num_threads(2)}


################################################
#### Load functions from the library ("Vw") ####
################################################

# First update the library and then load

#install.packages(Vw_library_path, repos = NULL, type = "source")
library(Vw)

finescale = TRUE # Whether a file containing data on locus-wise va_lost vs alpha should be saved
verbose = TRUE

########################
### Perform analyses ###
########################


Set_ID = ifelse(Sys.info()["nodename"]=="SCE-BIO-C06645"|Sys.info()["nodename"]=="sce-bio-c04553", "poo_seq_test_SCE-BIO-C06645_2025-07-04_12-53-13.067883", commandArgs(trailingOnly = TRUE)[1])
nsims = 1

for(sim in 1:nsims){

  sim_data = extract_slim_data(Set_ID = Set_ID,
                               sim = sim,
                               unzip = TRUE,
                               slim_output_path = slim_output_path, 
                               sim_param_path = sim_param_path,
                               extract_genomes_path = extract_genomes_path, 
                               extract_mut_path = extract_mut_path,
                               mutations_path = temp_files_path, 
                               c_matrix_path = temp_files_path, 
                               randomise = TRUE,
                               delete_temp_files = TRUE)
  
  c_genome = sim_data$c_genome  
  pbar0 = colMeans(c_genome)
  pbar2 = sim_data$pbar2
  diversity = pbar0*(1 - pbar0)/2
  n_cages = sim_data$sim_params$n_cages 
  n0_individuals = nrow(c_genome)/2
  
  list_alpha = sim_data$list_alpha
  
  # If the simulation was additive
  # "sim" in sim_data$sim_params does not contain information on dominance ("k")
  
  if(!grepl("k=", sim_data$sim_params$sim)){
    message("Computing the alphas...")
    list_alpha_new = list_alpha + 0.25*(1 - 2*pbar0)*list_alpha^2
  }else{ # If dominance related information is contained
    message("Computing the alphas with dominance...")
    
    k = as.numeric(unlist(strsplit(sim_data$sim_params$sim, "="))[2])
    
    # We want deleterious (beneficial) alleles to be recessive (dominance)
    # k must be positive (negative) when selectionCoeff is positive (negative). 
    
    k = k*sign(list_alpha)
    
    d = k*list_alpha
    qbar0 = 1 - pbar0
    list_alpha_new = list_alpha + d*(qbar0 - pbar0) + 0.25*(qbar0 - pbar0)*(list_alpha^2 + d*(2*list_alpha + d)) - 0.5*d*pbar0*((2*qbar0 - pbar0)*(2*list_alpha + qbar0*d) - pbar0*qbar0*d)
    
  }
  
  
  
  if(verbose){
    message("Calculating L...")
  }
  
  paternal<-seq(1, 2*n0_individuals, 2)
  maternal<-paternal+1
  
  c0<-(c_genome[paternal,]+c_genome[maternal,])/2
  
  rm("c_genome")
  
  L<-cov(c0)*(n0_individuals-1)/n0_individuals
  
  if(verbose){message("Calculating initial additive genetic and genic variances using old and new definitions of alphas...")}
  
  # Initial additive genetic variance
  
  vA_true = t(list_alpha)%*%L%*%list_alpha
  vA_true_new = t(list_alpha_new)%*%L%*%list_alpha_new
  
  # Initial additive genic variance
  
  va_true = sum(diversity*list_alpha^2)
  va_true_new = sum(diversity*list_alpha_new^2)
  
  if(verbose){message("Calculating average additive genic variance left at the end of the experiment using old and new definitions of alphas...")}
  
  
  # Average additive genic variance at the end of the experiment
  
  # Repeat list_alpha n_cages times over rows
  
  list_alpha_rep = t(matrix(list_alpha, nrow = length(list_alpha), ncol = n_cages))
  
  va_left = mean(rowSums(0.5*pbar2*(1-pbar2)*list_alpha_rep^2))
  
  list_alpha_rep_new = t(matrix(list_alpha_new, nrow = length(list_alpha_new), ncol = n_cages))
  
  va_left_new = mean(rowSums(0.5*pbar2*(1-pbar2)*list_alpha_rep_new^2))
  
  if(verbose){message("Saving data...")}
  
  ### Save file ###
  
  # Create a unique stamp for this analysis
  
  unique_stamp = as.character(paste(Sys.info()["nodename"], Sys.time()))
  unique_stamp = gsub(" ", "_", unique_stamp)
  unique_stamp = gsub(":", "-", unique_stamp)
  
  sim_params = sim_data$sim_params
  
  new_data = data.frame("vA_true" = vA_true, "vA_true_new" = vA_true_new, "va_true" = va_true, "va_true_new" = va_true_new, "va_left" = va_left, "va_left_new" = va_left_new, "time_stamp_va_left" = unique_stamp)
  new_data = cbind(sim_params, new_data)
  
  write.table(rbind(names(new_data), new_data), file = paste(output_path, "/", Set_ID, "_sim_", sim, "_new_va_calculation_", unique_stamp, ".csv", sep = ""),col.names = FALSE, row.names = FALSE, sep = ",")
  print(Set_ID)
  print(paste("Percentage additive genic variance left (average) = ", va_left/va_true*100, " (old) ", va_left_new/va_true_new*100, " (new) ", sep = ""))
  
  ### Optionally save fine-scale data ###
  # Save Set_ID, list_alpha, SNPs, locus-wise va_true and locus-wise va_left
  
  if(finescale){
    
    if(verbose){message("Saving finescale data...")}
    
    finescale_data = data.frame("Set_ID" = rep(Set_ID, length(list_alpha)),
                                 "sim" = rep(sim, length(list_alpha)),
                                 "list_alpha" = list_alpha,
                                 "list_alpha_new" = list_alpha_new,
                                 "pbar0" = pbar0,
                                 "pbar2_av" = colMeans(pbar2),
                                 "SNPs" = sim_data$SNPs,
                                 "locuswise_va_true" = diversity*list_alpha^2,
                                 "locuswise_va_left" = colMeans(0.5*pbar2*(1-pbar2)*list_alpha_rep^2),
                                 "locuswise_va_true_new" = diversity*list_alpha_new^2,
                                 "locuswise_va_left_new" = colMeans(0.5*pbar2*(1-pbar2)*list_alpha_rep_new^2),
                                 "vA_true" = rep(vA_true, length(list_alpha)),
                                 "vA_true_new" = rep(vA_true_new, length(list_alpha)),
                                 "va_true" = rep(va_true, length(list_alpha)),
                                 "va_true_new" = rep(va_true_new, length(list_alpha)))
    
    write.table(rbind(names(finescale_data), finescale_data), file = paste(output_path, "/", Set_ID, "_sim_", sim, "_finescale_va_", unique_stamp, ".csv", sep = ""),col.names = FALSE, row.names = FALSE, sep = ",")
    
  }
}