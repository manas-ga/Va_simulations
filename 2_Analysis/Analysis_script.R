
# Rscript /mnt/c/Users/msamant/Documents/GitHub/Va_simulations/2_Analysis/Analysis_script.R

# File Paths

if(Sys.info()["nodename"]!="SCE-BIO-C06645"){message("WARNING: Command line arguments required!!!")}

### On AC3 ###

if(Sys.info()["nodename"]%in%c("bigfoot", "bigshot", "bigbird", "bigyin", "biggar", "bigwig", "c1", "c2", "c3", "c4", "c5", "c6")){
  
  base_path = "/ceph/users/marun/Va_simulations/1_Simulations" 
  analysis_path = "/ceph/users/marun/Va_simulations/2_Analysis"
  file_storage_path = "/data/obbard/Va_simulations/analyses" # File storage path is the designated storage space on AC3 instead of the /home directory on qm
  
}

### On Vera ###

if(Sys.info()["nodename"]=="vera.bio.ed.ac.uk"){
  
  base_path = "/data/home/msamant/Manas/Va_simulations/Github/Va_simulations/1_Simulations" 
  analysis_path = "/data/home/msamant/Manas/Va_simulations/Github/Va_simulations/2_Analysis"  
  file_storage_path = base_path
  
}

### On Eddie ###

if(grepl("ecdf.ed.ac.uk", Sys.info()["nodename"])){
  base_path = "/home/msamant/Va_simulations/1_Simulations"
  analysis_path = "/home/msamant/Va_simulations/2_Analysis"
  file_storage_path = "/exports/eddie/scratch/msamant"
}

### On Manas's PC

if(Sys.info()["nodename"]=="SCE-BIO-C06645"){
  
  if(Sys.info()["sysname"]=="Linux"){   ## Local Wsl
    
    base_path = "/mnt/c/Users/msamant/Documents/GitHub/Va_simulations/1_Simulations"
    analysis_path = "/mnt/c/Users/msamant/Documents/GitHub/Va_simulations/2_Analysis" 
    file_storage_path = base_path
    
  }else{                                ## Local windows
    
    base_path = "C:/Users/msamant/Documents/GitHub/Va_simulations/1_Simulations" 
    analysis_path = "C:/Users/msamant/Documents/GitHub/Va_simulations/2_Analysis" 
    file_storage_path = base_path
  }
  
}

### On Jarrod's PC ###

if(Sys.info()["nodename"]=="sce-bio-c04553"){  
  base_path="~/Work/Va_simulations/1_Simulations"
  analysis_path = "~/Work/Va_simulations/2_Analysis"
  file_storage_path = base_path
  slim_path = "/Volumes/hadfield/Va_simulations/sim_files"
}else{
  slim_path = "/mnt/u/Datastore/CSCE/biology/groups/hadfield/Va_simulations/sim_files"
}


extract_genomes_path = file.path(analysis_path, "3_Extract_genomes.py")                                           ## Python script extracting mutations and genomes from the SLiM output file 
extract_mut_path = file.path(analysis_path, "2_Extract_mutations.py")                                             ## Python script extracting mutations from the SliM output file

slim_output_path = file.path(slim_path, "SLiM_outputs")
# Path to SLiM output              
sim_param_path = file.path(slim_path, "sim_params") 
# Path to SLiM simulation parameters

temp_files_path = file.path(file_storage_path, "temp_files")  
output_path = paste(file_storage_path, "/c_Output", sep = "")                                                     ## Path where .csv file(s) containing final data are to be stored 


# Load packages and functions

if(Sys.info()["nodename"]=="bigyin"){stop("Bigyin cannot run asreml-r. Use a different code.")}

library(MCMCglmm)
library(asreml)
library(Matrix)
library(pryr) ## For tracking memory usage using mem_used()
library(bigalgebra)
library(RhpcBLASctl)

# Control the number of BLAS threads if running on a cluster
if(Sys.info()["nodename"]!="SCE-BIO-C06645"|Sys.info()["nodename"]!="sce-bio-c04553"){blas_set_num_threads(10)}


#################################
#### Load Jarrod's functions ####
#################################

functions_only=TRUE ## Read only the functions

rmarkdown::render(file.path(analysis_path, "Vw.Rmd"))

##############################
### Load Manas's functions ###
##############################

rmarkdown::render(file.path(analysis_path, "Vw_sim_functions.Rmd"))

########################
### Perform analyses ###
########################

Set_ID = ifelse(Sys.info()["nodename"]=="SCE-BIO-C06645"|Sys.info()["nodename"]=="sce-bio-c04553", "test_unflipped_SCE-BIO-C06645_2024-09-11_13_24_45.087551", commandArgs(trailingOnly = TRUE)[1])
sim = ifelse(Sys.info()["nodename"]=="SCE-BIO-C06645"|Sys.info()["nodename"]=="sce-bio-c04553", 1, as.numeric(commandArgs(trailingOnly = TRUE)[1]))

analysed_data = analyse_sim(Set_ID = Set_ID,                            # The unique ID of the set of simulations that are controlled by a single R script
                            sim = sim,                                    # Each set can have multiple sims, but - on the cluster sim must always 1
                            unzip = TRUE,                              # Should the SLiM output file be unzipped, read, and then zipped back?
                            slim_output_path = slim_output_path,        # The directory where the SLiM outputs (for parents and experimental replicates) are stored (as .txt files)
                            sim_param_path = sim_param_path,            # The path to the directory where the .csv file containing simulation parameters is stored
                            extract_genomes_path = extract_genomes_path,# The path to the python script that extracts genomes and mutations from SLim outputs
                            extract_mut_path = extract_mut_path,        # The path to the python script that extracts mutations from SLim outputs
                            mutations_path = temp_files_path,           # The directory where extracted mutations are to be stored (temp files)
                            c_matrix_path = temp_files_path,            # The directory where extracted genomes are to be stored (temp files)
                            output_path = output_path,                  # The path where the final data file is to be stored
                            n_sample=NULL,                              # Number of individuals sampled from the parents' generation (useful if n_ind_exp is large)
                            randomise = TRUE,                           # Optionally the reference allele can be randomised
                            delete_temp_files = TRUE,
                            proj = "BLoM",                              # projection type for allele frequencies: "LoM", "BLoM", "L" or "N"
                            LDdelta = FALSE,                            # Should L or diag(L) be considered while modelling distribution of alphas
                            pa = 1,
                            Vs = "LoNL",                                # "L" or "LoNL"
                            method="REML",                              # Can be "REML" or "MCMC"
                            pdelta = NA,                                # If NA pdelta is estimated using optim()
                            bdelta = c(NA, NA),                         # If c(NA,NA) both bedelta intercept and slope are estimated
                            AtleastOneRecomb=FALSE,
                            verbose = TRUE)