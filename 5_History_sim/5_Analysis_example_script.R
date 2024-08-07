
# source("/mnt/c/Users/msamant/Documents/GitHub/Va_simulations/5_History_sim/5_Analysis_example_script.R")

##############################################################
########## Read, process, and analyse SLiM outputs ###########
##############################################################

rm(list = ls())

########################################################################
########### paths of various scripts and functions #####################
########################################################################

if(Sys.info()["nodename"]!="SCE-BIO-C06645"){message("WARNING: Command line arguments required!!!")}

### Base path and path to Vw.Rmd (file containing Jarrod's functions) (depending on the system) ###

### On AC3 ###

if(Sys.info()["nodename"]%in%c("bigfoot", "bigshot", "bigbird", "bigyin", "biggar", "bigwig", "c1", "c2", "c3", "c4", "c5", "c6")){
  
  base_path = "/ceph/users/marun/Va_simulations/5_History_sim" # Path to all the scripts except "Vw.Rmd"
  Vw_path = "/ceph/users/marun/Va_simulations/6_Code_Test/Vw.Rmd"
  file_storage_path = "/data/obbard/Va_simulations/analyses" # File storage path is the designated storage space on AC3 instead of the /home directory on qm
  
}

### On Vera ###

if(Sys.info()["nodename"]=="vera.bio.ed.ac.uk"){
  
  base_path = "/data/home/msamant/Manas/Va_simulations/Github/Va_simulations/5_History_sim" ## ON VERA
  Vw_path = "/data/home/msamant/Manas/Va_simulations/Github/Va_simulations/6_Code_Test/Vw.Rmd"  ### Jarrod's functions and other code is stored here
  file_storage_path = base_path
  
}

### On Eddie ###

if(grepl("ecdf.ed.ac.uk", Sys.info()["nodename"])){
  base_path = "/home/msamant/Va_simulations/5_History_sim"
  Vw_path = "/home/msamant/Va_simulations/6_Code_Test/Vw.Rmd"
  file_storage_path = "/exports/eddie/scratch/msamant"
}

### On Manas's PC

if(Sys.info()["nodename"]=="SCE-BIO-C06645"){
  
  if(Sys.info()["sysname"]=="Linux"){
    
    base_path = "/mnt/c/Users/msamant/Documents/GitHub/Va_simulations/5_History_sim" ## Local Wsl
    Vw_path = "/mnt/c/Users/msamant/Documents/GitHub/Va_simulations/6_Code_test/Vw.Rmd" ### Jarrod's functions and other code is stored here
    file_storage_path = base_path
    
  }else{
    
    base_path = "C:/Users/msamant/Documents/GitHub/Va_simulations/5_History_sim" ## Local windows
    Vw_path = "C:/Users/msamant/Documents/GitHub/Va_simulations/6_Code_test/Vw.Rmd" ### Jarrod's functions and other code is stored here
    file_storage_path = base_path
  }
  
}

### On Jarrod's PC ###

if(Sys.info()["nodename"]=="sce-bio-c04553"){  
  base_path="~/Work/Va_simulations//5_History_sim"
  Vw_path = "~/Work/Va_simulations/6_Code_test/Vw.Rmd" ### Jarrod's functions and other code is stored here
  file_storage_path = base_path
}

# Paths to various scripts that are used for running the simulations and extracting information from SLiM outputs

extract_genomes_path = file.path(base_path, "3_Extract_genomes.py")                                           ## Python script extracting mutations and genomes from the SLiM output file 
extract_mut_path = file.path(base_path, "2_Extract_mutations.py")                                             ## Python script extracting mutations from the SliM output file

# Paths to various directories !!!!!! NEEDS TO BE SPECIFIED !!!!!!

slim_output_path = "/mnt/u/Datastore/CSCE/biology/groups/hadfield/Va_simulations/sim_files/SLiM_outputs"                                 ## Path where SLiM and msprime output files are stored
slim_param_path = "/mnt/u/Datastore/CSCE/biology/groups/hadfield/Va_simulations/sim_files/sim_params"
temp_files_path = paste(file_storage_path, "/temp_files", sep = "")  

####################################
######### Packages #################
####################################

if(Sys.info()["nodename"]=="bigyin"){stop("Bigyin cannot run asreml-r. Use a different code.")}

library(MCMCglmm)
library(asreml)
library(Matrix)
library(rmutil)
library(pryr) ## For tracking memory usage using mem_used()
library(bigalgebra)
library(RhpcBLASctl)

  
#################################
#### Load Jarrod's functions ####
#################################

functions_only=TRUE ## Read only the functions

rmarkdown::render(file.path(Vw_path))

##############################
### Load Manas's functions ###
##############################

rmarkdown::render(file.path(base_path, "Vw_sim_functions.Rmd"))

#### Enter the Set_ID of the simulations to be analysed ###

Set_ID = "bigfoot_2024-08-01_13:03:32.727659_9.64489795918367e-07_1.4_1.4_1000_10_3_estimate_0"

### Extract sim data ###

sim_data = extract_slim_data(Set_ID = Set_ID,
                             sim = 1,
                             slim_output_path = slim_output_path, 
                             sim_param_path = slim_param_path,
                             extract_genomes_path = extract_genomes_path, 
                             extract_mut_path = extract_mut_path,
                             mutations_path = temp_files_path, 
                             c_matrix_path = temp_files_path, 
                             randomise = TRUE)
### Analyse ###

parents_info = analyse_parents(c0 = sim_data$c0,  
                               list_alpha = sim_data$list_alpha,             
                               LDdelta=FALSE,         
                               SNPs = sim_data$SNPs,                   
                               RecombRate = sim_data$sim_params$r_expt,             
                               HapLength = sim_data$sim_params$sequence_length,              
                               AtleastOneRecomb=FALSE)

### Fit model ###

# Analysis parameters
proj="BLoM" # projection type for allele frequencies: "LoM", "BLoM", "L" or "N"
LDdelta = FALSE
pa = 1
Vs = "LoNL" # "L" or "LoNL"
method="REML" # Can be "REML" or "MCMC"
randomise = TRUE # Should the reference allele be randomised for analysis?
bigalgebra = FALSE # Should bigalgebra be used for eigendecomposition?

# How is pdelta to be estimated? 
# Can be "optim" (using the function optim()), or "fixed" or "manual"(estimated by manually scanning a range of pdelta values)

pdelta_method = "optim" # "optim" or "manual" or "fixed" or "no_analysis". If this is "no_analysis", the estimate of Vw is not calculated, but the rest of the code still runs.

if(pdelta_method=="fixed"){
  pdelta = 0 # Can be specified to any value
}

if(pdelta_method=="optim"){
  pdelta = NA # This triggers the use of optim() inside the function Vw_model()
}

if(pdelta_method=="manual"){
  
  nseq<-20 # The number of times pdelta is to be varied 
  pdelta_l = -1 # Lower limit of pdelta
  pdelta_u = 1 # Upper limit of pdelta
  pdelta<-seq(pdelta_l, pdelta_u, length=nseq)
  
}

# How should bdelta[1] (intercept) and bdelta[2] (slope of (p-q)) be estimated

bdelta_method = ifelse(Sys.info()["nodename"]=="SCE-BIO-C06645", "estimate", commandArgs(trailingOnly = TRUE)[7])  # Can be "fixed" or "estimate"

if(bdelta_method=="estimate"){
  bdelta = c(NA, NA)
}else{
  bdelta = c(0, 0) # This only estimates the intercept while keeping the slope fixed at 0
}


m1<-Vw_model(C0 = sim_data$c0,          # parental genotypes (rows individuals, columns loci, coded as 0, 1/2 or 1) 
             nR = parents_info$nR,          # matrix of non-recombinant probabilities between loci
             pbar1 = sim_data$pbar1,       # vector of allele frequencies at time-point 1
             ngen1=1,     # number of generations between parents and time-point 1
             pbar2 = sim_data$pbar2,       # vector of allele frequencies at time-point 2
             ngen2 = (sim_data$sim_params$ngen_expt) + 1,       # number of generations between parents and time-point 2
             nind = sim_data$sim_params$n_ind_exp,        # population size in each replicate
             proj=proj, # projection type for allele frequencies: "LoM", "BLoM", "L" or "N"
             LDdelta = LDdelta,
             pa = pa,
             pdelta = pdelta,
             bdelta = bdelta,
             Vs = Vs,
             method = method,
             L = parents_info$L,    # list with elements UL and DL
             svdL = NULL,    # list with elements UL and DL
             bigalgebra = bigalgebra, 
             tol = sqrt(.Machine$double.eps))


vA_est = m1$Vw_est # Store the estimate from the model with the highest log likelihood
pdelta_est = m1$pdelta # The sample() functions ensures that only one value is selected in case there are multiple points with the highest LL
pdelta_var_est = m1$pdelta_var
bdelta_intercept_est = m1$bdelta[1]
bdelta_slope_est = m1$bdelta[2]
bdelta_var_est = paste(m1$bdelta_var[1,1], m1$bdelta_var[2,2], m1$bdelta_var[1,2], sep = "_")
sigma2delta_est = summary(m1$model)$varcomp[1,1]



