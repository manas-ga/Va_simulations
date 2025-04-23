# This script calculates V_A based on Buffalo and Coop's (2019) method in three different ways
# 1. del_P for all segregating sites and average LD using all segregating sites
# 2. del_P for neutral sites and average LD using all segregating sites
# 3. del_P for neutral sites and average LD using the LD between selected and neutral sites only

rm(list = ls())
# Rscript /mnt/c/Users/msamant/Documents/GitHub/Va_simulations/2_Analysis/bc_va_script.R
# source("/mnt/c/Users/msamant/Documents/GitHub/Va_simulations/2_Analysis/bc_va_script.R")


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
if(Sys.info()["nodename"]!="SCE-BIO-C06645"|Sys.info()["nodename"]!="sce-bio-c04553"){blas_set_num_threads(10)}


################################################
#### Load functions from the library ("Vw") ####
################################################

# First update the library and then load

#install.packages(Vw_library_path, repos = NULL, type = "source")
library(Vw)

verbose = TRUE

########################
### Perform analyses ###
########################


Set_ID = ifelse(Sys.info()["nodename"]=="SCE-BIO-C06645"|Sys.info()["nodename"]=="sce-bio-c04553", "bc_test_SCE-BIO-C06645_2025-04-22_16-02-10.727186", commandArgs(trailingOnly = TRUE)[1])
nsims = 1

for(sim in 1:nsims){
  
  sim_data = extract_slim_data(Set_ID = Set_ID,
                               sim = sim,
                               ngen2_optional = 2,
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
  n0_individuals = nrow(c_genome)/2
  pbar0 = colMeans(c_genome)
  diversity = pbar0*(1 - pbar0)/2
  list_alpha = sim_data$list_alpha
  list_alpha_new = list_alpha + 0.25*(1 - 2*pbar0)*list_alpha^2
  
  if(verbose){
    message("Calculating L...")
  }
  
  paternal<-seq(1, 2*n0_individuals, 2)
  maternal<-paternal+1
  
  c0<-(c_genome[paternal,]+c_genome[maternal,])/2
  
  rm("c_genome")
  
  L<-cov(c0)*(n0_individuals-1)/n0_individuals
  
  rm("c0")
  
  
  message("Calculating the true levels of V_A...")
  
  # Initial additive genetic variance
  
  vA_true = t(list_alpha)%*%L%*%list_alpha
  vA_true_new = t(list_alpha_new)%*%L%*%list_alpha_new
  
  # Initial additive genic variance
  
  va_true = sum(diversity*list_alpha^2)
  va_true_new = sum(diversity*list_alpha_new^2)
  
  
  
  if(verbose){
    message("Calculating the matrix of non-recombinant fractions...")
  }
  
  nR = form_nR(SNPs = sim_data$SNPs, RecombRate = sim_data$sim_params$r_expt, HapLength = sim_data$sim_params$sequence_length, AtleastOneRecomb = FALSE)

  ### Calculate Vw from Buffalo and Coop's method ###
  
  if(sim_data$ngen2-sim_data$ngen1==1){
    
    message("Calculating Vw using Buffalo and Coop's (2019) method (approach 1) ...")
    
    ### Approach 1: del_P for all segregating sites and average LD using all segregating sites
    
    BC_fit_1 = est_Va_bc(pbar1 = sim_data$pbar1,
                         pbar2 = sim_data$pbar2,
                         L = L,
                         nR = nR)
    
    ### Approach 2: del_P for neutral sites and average LD using all segregating sites
    message("Calculating Vw using Buffalo and Coop's (2019) method (approach 2) ...")
    # Identify selected loci
    
    selected = which(sim_data$list_alpha!=0)
    
    BC_fit_2 = est_Va_bc(pbar1 = sim_data$pbar1[,-selected],
                         pbar2 = sim_data$pbar2[,-selected],
                         L = L,
                         nR = nR)
    
    
    # Approach 3: del_P for neutral sites and average LD using the LD between selected and neutral sites only
    message("Calculating Vw using Buffalo and Coop's (2019) method (approach 3) ...")
    BC_fit_3 = est_Va_bc(pbar1 = sim_data$pbar1[,-selected],
                         pbar2 = sim_data$pbar2[,-selected],
                         L = L,
                         nR = nR,
                         selected = selected)
    
  }
  
  if(verbose){message("Saving data...")}
  
  #################
  ### Save data ###
  #################
  
  # Create a unique stamp for this analysis
  
  unique_stamp = as.character(paste(Sys.info()["nodename"], Sys.time()))
  unique_stamp = gsub(" ", "_", unique_stamp)
  unique_stamp = gsub(":", "-", unique_stamp)
  
  sim_params = sim_data$sim_params
  
  new_data = data.frame("vA_true" = vA_true, "vA_true_new" = vA_true_new, "va_true" = va_true, "va_true_new" = va_true_new,
                        "vA_BC_1" = if(sim_data$ngen2-sim_data$ngen1==1){BC_fit_1$vA_BC}else{NA}, "Ne_BC_1" = if(sim_data$ngen2-sim_data$ngen1==1){BC_fit_1$Ne_BC}else{NA},
                        "vA_BC_2" = if(sim_data$ngen2-sim_data$ngen1==1){BC_fit_2$vA_BC}else{NA}, "Ne_BC_2" = if(sim_data$ngen2-sim_data$ngen1==1){BC_fit_2$Ne_BC}else{NA},
                        "vA_BC_3" = if(sim_data$ngen2-sim_data$ngen1==1){BC_fit_3$vA_BC}else{NA}, "Ne_BC_3" = if(sim_data$ngen2-sim_data$ngen1==1){BC_fit_3$Ne_BC}else{NA},
                        "time_stamp_va_left" = unique_stamp)
  new_data = cbind(sim_params, new_data)
  
  write.table(rbind(names(new_data), new_data), file = paste(output_path, "/", Set_ID, "_sim_", sim, "_bc_va_corrected_", unique_stamp, ".csv", sep = ""),col.names = FALSE, row.names = FALSE, sep = ",")
  message(paste("V_A estimated using the Buffalo & Coop (2019) method for", Set_ID, "sim", sim))
  
}
