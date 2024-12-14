rm(list = ls())
# Rscript /mnt/c/Users/msamant/Documents/GitHub/Va_simulations/2_Analysis/Analysis_script.R
# source("/mnt/c/Users/msamant/Documents/GitHub/Va_simulations/2_Analysis/Analysis_script.R")

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
if(Sys.info()["nodename"]!="SCE-BIO-C06645"|Sys.info()["nodename"]!="sce-bio-c04553"){blas_set_num_threads(15)}


################################################
#### Load functions from the library ("Vw") ####
################################################

# First update the library and then load

#install.packages(Vw_library_path, repos = NULL, type = "source")
library(Vw)

########################
### Perform analyses ###
########################


Set_ID = ifelse(Sys.info()["nodename"]=="SCE-BIO-C06645"|Sys.info()["nodename"]=="sce-bio-c04553", "noburnin_test_SCE-BIO-C06645_2024-12-13_23-22-02.930194", commandArgs(trailingOnly = TRUE)[1])
nsims = ifelse(Sys.info()["nodename"]=="SCE-BIO-C06645"|Sys.info()["nodename"]=="sce-bio-c04553", 1, as.numeric(commandArgs(trailingOnly = TRUE)[2]))

for(sim in 1:nsims){
  message(paste("Analysing simulation", sim, "of set", Set_ID, "..."))
  analysed_data = analyse_sim(Set_ID = Set_ID,                            # The unique ID of the set of simulations that are controlled by a single R script
                              sim = sim,                                  # Each set can have multiple sims, but - on the cluster sim must always 1
                              unzip = TRUE,                               # Should the SLiM output file be unzipped, read, and then zipped back?
                              slim_output_path = slim_output_path,        # The directory where the SLiM outputs (for parents and experimental replicates) are stored (as .txt files)
                              sim_param_path = sim_param_path,            # The path to the directory where the .csv file containing simulation parameters is stored
                              extract_genomes_path = extract_genomes_path,# The path to the python script that extracts genomes and mutations from SLim outputs
                              extract_mut_path = extract_mut_path,        # The path to the python script that extracts mutations from SLim outputs
                              mutations_path = temp_files_path,           # The directory where extracted mutations are to be stored (temp files)
                              c_matrix_path = temp_files_path,            # The directory where extracted genomes are to be stored (temp files)
                              output_path = output_path,                  # The path where the final data file is to be stored
                              randomise = TRUE,                           # Optionally the reference allele can be randomised
                              delete_temp_files = TRUE,
                              proj = "BLoM",                              # projection type for allele frequencies: "LoM", "BLoM", "L" or "N"
                              LDalpha = FALSE,                            # Should L or diag(L) be considered while modelling distribution of alphas
                              pa = 1,
                              Vs = "LoNL",                                # "L" or "LoNL"
                              method="REML",                              # Can be "REML" or "MCMC"
                              palpha = NA,                                # If NA pdelta is estimated using optim()
                              balpha = c(NA, NA),                         # If c(NA,NA) both bedelta intercept and slope are estimated
                              AtleastOneRecomb=FALSE,
                              Ne = c(1000, 1000),
                              predict_Ne =  TRUE,                         # If true, this overwrites the Ne supplied above by Ne = c(nind_expt, predict_Ne(nind_expt, Ve_w_expt))
                              verbose = TRUE)
  
  message("Final output:")
  print(analysed_data)

}