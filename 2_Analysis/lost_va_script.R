##################
### File Paths ###
##################

if(Sys.info()["nodename"]!="SCE-BIO-C06645"){message("WARNING: Command line arguments required!!!")}

### On AC3 ###

if(Sys.info()["nodename"]%in%c("bigfoot", "bigshot", "bigbird", "bigyin", "biggar", "bigwig", "c1", "c2", "c3", "c4", "c5", "c6")){
  analysis_path = "/ceph/users/marun/Va_simulations/2_Analysis"
  file_storage_path = "/data/obbard/Va_simulations/analyses" # File storage path is the designated storage space on AC3 instead of the /home directory on qm
  Vw_library_path = "/ceph/users/marun/Va_simulations/Vw"
}

### On Vera ###

if(Sys.info()["nodename"]=="vera.bio.ed.ac.uk"){
  analysis_path = "/data/home/msamant/Manas/Va_simulations/Github/Va_simulations/2_Analysis"  
  file_storage_path = "/data/home/msamant/Manas/Va_simulations/Github/Va_simulations/1_Simulations"
  Vw_library_path = "/data/home/msamant/Manas/Va_simulations/Github/Va_simulations/Vw"
  
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

if(Sys.info()["nodename"]%in%c("bigfoot", "bigshot", "bigbird", "bigyin", "biggar", "bigwig", "c1", "c2", "c3", "c4", "c5", "c6")){
  slim_output_path = file.path(file_storage_path, "b_Interim_files/SLiM_outputs")
}

##########################################
### Path to SLiM simulation parameters ###
##########################################

if(Sys.info()["nodename"]=="SCE-BIO-C06645"|Sys.info()["nodename"]=="sce-bio-c04553"){
  sim_param_path = if(test){file.path(file_storage_path, "c_Output")}else{file.path(slim_path, "sim_params")}   
}

if(Sys.info()["nodename"]%in%c("bigfoot", "bigshot", "bigbird", "bigyin", "biggar", "bigwig", "c1", "c2", "c3", "c4", "c5", "c6")){
  sim_param_path = file.path(file_storage_path, "c_Output")
}

# Load packages and functions

if(Sys.info()["nodename"]=="bigyin"){stop("Bigyin cannot run asreml-r. Use a different node.")}

library(MCMCglmm)
library(asreml)
library(Matrix)
library(pryr) ## For tracking memory usage using mem_used()
library(RhpcBLASctl)

# Control the number of BLAS threads if running on a cluster
if(Sys.info()["nodename"]!="SCE-BIO-C06645"|Sys.info()["nodename"]!="sce-bio-c04553"){blas_set_num_threads(22)}


################################################
#### Load functions from the library ("Vw") ####
################################################

# First update the library and then load

#install.packages(Vw_library_path, repos = NULL, type = "source")
library(Vw)


########################
### Perform analyses ###
########################


Set_ID = ifelse(Sys.info()["nodename"]=="SCE-BIO-C06645"|Sys.info()["nodename"]=="sce-bio-c04553", "Set_5_c3_2024-12-17_02-56-59.533089_4.65656565656566e-09_5_2_1000_10_4_0_1", commandArgs(trailingOnly = TRUE)[1])
nsims = ifelse(Sys.info()["nodename"]=="SCE-BIO-C06645"|Sys.info()["nodename"]=="sce-bio-c04553", 1, as.numeric(commandArgs(trailingOnly = TRUE)[2]))


for(sim in 1:nsims){
  
  sim_data = extract_slim_data(Set_ID = Set_ID,
                               sim = sim,
                               unzip = unzip,
                               slim_output_path = slim_output_path, 
                               sim_param_path = sim_param_path,
                               extract_genomes_path = extract_genomes_path, 
                               extract_mut_path = extract_mut_path,
                               mutations_path = mutations_path, 
                               c_matrix_path = c_matrix_path, 
                               randomise = TRUE,
                               delete_temp_files = delete_temp_files)
  
  c_genome = sim_data$c_genome  
  list_alpha = sim_data$list_alpha
  pbar2 = sim_data$pbar2
  n_cages = sim_data$sim_params$n_cages 
  
  # Initial additive genic variance
  
  pbar0 = colMeans(c_genome) 
  diversity = pbar0*(1 - pbar0)/2
  va_true = sum(diversity*list_alpha^2)
  
  # Average additive genic variance at the end of the experiment
  
  # Repeat list_alpha n_cages times over rows
  
  list_alpha = t(matrix(list_alpha, nrow = length(list_alpha), ncol = n_cages))
  
  va_left = mean(rowSums(0.5*pbar2*(1-pbar2)*list_alpha))
  
  
}