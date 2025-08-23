rm(list = ls())
# Jarrod's PC ------ "sce-bio-c04553"
# Rscript ~/Work/Va_simulations/2_Analysis/Qpoolseq_diag_script.R
# source("~/Work/Va_simulations/2_Analysis/Qpoolseq_diag_script.R")

# Manas's PC ------ "SCE-BIO-C06645"
# Rscript /mnt/c/Users/msamant/Documents/GitHub/Va_simulations/2_Analysis/Qpoolseq_diag_script.R
# source("/mnt/c/Users/msamant/Documents/GitHub/Va_simulations/2_Analysis/Qpoolseq_diag_script.R")

##################
### File Paths ###
##################

if(Sys.info()["nodename"]!="SCE-BIO-C06645"){message("WARNING: Command line arguments required!!!")}

### On AC3 ###

if(Sys.info()["nodename"]%in%c("bigfoot", "bigshot", "bigbird", "bigyin", "biggar", "bigwig", "c1", "c2", "c3", "c4", "c5", "c6", "ac3-n1", "ac3-n2", "ac3-n3", "ac3-n4", "ac3-n5", "ac3-n6")){
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
  file_storage_path = "/Volumes/hadfield/Va_simulations/sim_files"
  Vw_library_path = "~/Work/Va_simulations/Vw"
}

####################
### Script paths ###
####################

extract_genomes_path = file.path(analysis_path, "3_Extract_genomes.py")                                           ## Python script extracting mutations and genomes from the SLiM output file 
extract_mut_path = file.path(analysis_path, "2_Extract_mutations.py")                                             ## Python script extracting mutations from the SliM output file

#########################
### Output file paths ###
#########################

temp_files_path = file.path(file_storage_path, "b_Interim_files")  
output_path = file.path(file_storage_path, "c_Output")                                                     ## Path where .csv file(s) containing final data are to be stored 

# If working locally on Manas' or Jarrod's PC, send the outputs to local destinations (not Datastore)

if(Sys.info()["nodename"]=="SCE-BIO-C06645"|Sys.info()["nodename"]=="sce-bio-c04553"){
  temp_files_path = file.path(analysis_path, "Output")  
  output_path = file.path(analysis_path, "Output")   
}

# On Eddie use a different output path (since files in /scracth are deleted after a months)

if(grepl("ecdf.ed.ac.uk", Sys.info()["nodename"])){
  output_path = paste("~", "/c_Output", sep = "")  
}

###########################
### Path to SLiM output ###
###########################

test = TRUE # If TRUE, this allows Manas to test stuff offline using sims stored on his PC (as opposed to using sims on Datastore)

if(Sys.info()["nodename"]=="SCE-BIO-C06645"){
  slim_output_path = if(test){file.path(file_storage_path, "b_Interim_files/SLiM_outputs")}else{file.path(slim_path, "SLiM_outputs")}
}

if(Sys.info()["nodename"]=="sce-bio-c04553"){
  slim_output_path = file.path(file_storage_path, "b_Interim_files/SLiM_outputs")
}

if(Sys.info()["nodename"]%in%c("bigfoot", "bigshot", "bigbird", "bigyin", "biggar", "bigwig", "c1", "c2", "c3", "c4", "c5", "c6", "ac3-n1", "ac3-n2", "ac3-n3", "ac3-n4", "ac3-n5", "ac3-n6")){
  slim_output_path = file.path(file_storage_path, "b_Interim_files/SLiM_outputs")
}

if(grepl("ecdf.ed.ac.uk", Sys.info()["nodename"])){
  slim_output_path = file.path(file_storage_path, "b_Interim_files/SLiM_outputs")
}



##########################################
### Path to SLiM simulation parameters ###
##########################################

if(Sys.info()["nodename"]=="SCE-BIO-C06645"){
  sim_param_path = if(test){file.path(file_storage_path, "c_Output")}else{file.path(slim_path, "sim_params")}   
}

if(Sys.info()["nodename"]=="sce-bio-c04553"){
  sim_param_path = file.path(file_storage_path, "sim_params")
}

if(Sys.info()["nodename"]%in%c("bigfoot", "bigshot", "bigbird", "bigyin", "biggar", "bigwig", "c1", "c2", "c3", "c4", "c5", "c6", "ac3-n1", "ac3-n2", "ac3-n3", "ac3-n4", "ac3-n5", "ac3-n6")){
  sim_param_path = file.path(file_storage_path, "c_Output")
}

if(grepl("ecdf.ed.ac.uk", Sys.info()["nodename"])){
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
if(Sys.info()["nodename"]!="SCE-BIO-C06645"|Sys.info()["nodename"]!="sce-bio-c04553"){blas_set_num_threads(2)}


################################################
#### Load functions from the library ("Vw") ####
################################################

# First update the library and then load

#install.packages(Vw_library_path, repos = NULL, type = "source")
library(Vw)

### Set_ID ###

Set_ID = ifelse(Sys.info()["nodename"]=="SCE-BIO-C06645"|Sys.info()["nodename"]=="sce-bio-c04553", "TEST_SCE-BIO-C06645_2025-08-21_16-01-17.50597", commandArgs(trailingOnly = TRUE)[1])
sim = 1

### Pool-seq paramterers ###

pool_seq = TRUE

if(pool_seq){
  asreml.options(workspace="4gb") # only for poolseq
  read_length = 37
  coverage = 50
  V_logmean = 0
}else{
  read_length = NULL
  coverage = NULL
  V_logmean = NULL
}


### Analysis params ###

ngen2_optional = NULL                      # Allows del_P to be calculated between ngen1 and manually specified ngen2 (which can be different from the last generation)
unzip = TRUE                               # Should the SLiM output file be unzipped, read, and then zipped back?
slim_output_path = slim_output_path        # The directory where the SLiM outputs (for parents and experimental replicates) are stored (as .txt files)
sim_param_path = sim_param_path            # The path to the directory where the .csv file containing simulation parameters is stored
extract_genomes_path = extract_genomes_path# The path to the python script that extracts genomes and mutations from SLim outputs
extract_mut_path = extract_mut_path        # The path to the python script that extracts mutations from SLim outputs
mutations_path = temp_files_path           # The directory where extracted mutations are to be stored (temp files)
c_matrix_path = temp_files_path            # The directory where extracted genomes are to be stored (temp files)
output_path = output_path                  # The path where the final data file is to be stored
randomise = TRUE                           # Optionally the reference allele can be randomised
delete_temp_files = TRUE
pool_seq = pool_seq                        # Should the function simulate_pool_seq be used to sample allele frequencies in the experiment?
read_length = read_length
coverage = coverage
V_logmean = V_logmean
proj = "BLoM"                              # projection type for allele frequencies: "LoM", "BLoM", "L" or "N"
LDalpha = FALSE                            # Should L or diag(L) be considered while modelling distribution of alphas
pa = 1
Vs = "LoNL"                                # "L" or "LoNL"
method="REML"                              # Can be "REML" or "MCMC"
palpha = 0                                # If NA pdelta is estimated using optim()
balpha = c(0, NA)                          # If c(NA,NA) both bedelta intercept and slope are estimated
AtleastOneRecomb=FALSE
NE = c(1000, 1000)
Ne=NE
predict_NE =  TRUE                         # If true, this overwrites the Ne supplied above by Ne = c(nind_expt, predict_Ne(nind_expt, Ve_w_expt))
verbose = TRUE
all.gp = FALSE


### Checks

if(pool_seq){
  if(is.null(read_length)|is.null(coverage)|is.null(V_logmean)){stop("If poos_seq is true, read_length, coverage, and V_logmean must be provided")}
}

### Extract sim data ###

if(verbose){message("Extracting simulation data...")}

sim_data = extract_slim_data(Set_ID = Set_ID,
                             sim = sim,
                             ngen2_optional = ngen2_optional,
                             unzip = unzip,
                             slim_output_path = slim_output_path, 
                             sim_param_path = sim_param_path,
                             extract_genomes_path = extract_genomes_path, 
                             extract_mut_path = extract_mut_path,
                             mutations_path = mutations_path, 
                             c_matrix_path = c_matrix_path, 
                             randomise = randomise,
                             delete_temp_files = delete_temp_files,
                             pool_seq = pool_seq,
                             read_length = read_length,
                             coverage = coverage,
                             V_logmean = V_logmean)
### Analyse parents ###

if(verbose){message("Analysing the parents' generation...")}

parents_info = analyse_parents(c_genome = sim_data$c_genome,  
                               list_alpha = sim_data$list_alpha,     
                               compute_svdL=TRUE,        
                               LDalpha=LDalpha,   
                               SNPs = sim_data$SNPs,                   
                               RecombRate = sim_data$sim_params$r_expt,             
                               HapLength = sim_data$sim_params$sequence_length,              
                               AtleastOneRecomb=AtleastOneRecomb)

### Predict Ne in the experiment ###

if(predict_NE){
  NE<- c(sim_data$sim_params$n_ind_exp, predict_NE(n=sim_data$sim_params$n_ind_exp, Ve=sim_data$sim_params$Ve_w_expt))
}


### Fit model with diagonal Q ###

if(verbose){message("Performing analyses with diagonal Q...")}

m1_diag<-Vw_model(c_genome = NULL,          
             nR = parents_info$nR,
             pbar0 = sim_data$pbar0,                   
             pbar1 = sim_data$pbar1,
             coverage1 = sim_data$coverage1,
             ngen1=sim_data$ngen1,     
             pbar2 = sim_data$pbar2,  
             coverage2 = sim_data$coverage2,
             ngen2 = sim_data$ngen2,       
             nind = sim_data$sim_params$n_ind_exp,        
             proj=proj,
             LDalpha = LDalpha,
             pa = pa,
             palpha = palpha,
             balpha = balpha,
             Vs = Vs,
             method = method,
             L = parents_info$L,
             Ltilde = if(all.gp){parents_info$L}else{parents_info$Ltilde},      
             svdL = parents_info$svdL,           # list with elements UL and DL
             Ne = Ne,
             NE = NE,
             tol = sqrt(.Machine$double.eps))


vA_est_diag = m1_diag$Vw_est 
palpha_est_diag = m1_diag$palpha 
palpha_var_est_diag = m1_diag$palpha_var
balpha_intercept_est_diag = m1_diag$balpha[1]
balpha_slope_est_diag = m1_diag$balpha[2]
balpha_var_est_diag = paste(m1_diag$balpha_var[1,1], m1_diag$balpha_var[2,2], m1_diag$balpha_var[1,2], sep = "_")
sigma2alpha_est_diag = summary(m1_diag$model)$varcomp['vm(locus, SC, singG = "PSD")', 'component'] # change to summary(m1_diag$model)$varcomp['vm(locus, SC, singG = "PSD")', 'component']
Residual_var_diag = summary(m1_diag$model)$varcomp['units!R', 'component']

### Fit model with complete Q ###

if(verbose){message("Performing analyses with complete Q...")}

m1_comp<-Vw_model(c_genome = NULL,          
             nR = parents_info$nR,
             pbar0 = sim_data$pbar0,                   
             pbar1 = sim_data$pbar1,
             coverage1 = sim_data$coverage1,
             ngen1=sim_data$ngen1,     
             pbar2 = sim_data$pbar2,  
             coverage2 = sim_data$coverage2,
             ngen2 = sim_data$ngen2,       
             nind = sim_data$sim_params$n_ind_exp,        
             proj=proj,
             LDalpha = LDalpha,
             pa = pa,
             palpha = palpha,
             balpha = balpha,
             Vs = Vs,
             method = method,
             L = parents_info$L,
             Ltilde = if(all.gp){parents_info$L}else{parents_info$Ltilde},      
             svdL = parents_info$svdL,           # list with elements UL and DL
             Ne = Ne,
             NE = NE,
             diag.projQ=FALSE,
             tol = sqrt(.Machine$double.eps))


vA_est_comp = m1_comp$Vw_est 
palpha_est_comp = m1_comp$palpha 
palpha_var_est_comp = m1_comp$palpha_var
balpha_intercept_est_comp = m1_comp$balpha[1]
balpha_slope_est_comp = m1_comp$balpha[2]
balpha_var_est_comp = paste(m1_comp$balpha_var[1,1], m1_comp$balpha_var[2,2], m1_comp$balpha_var[1,2], sep = "_")
sigma2alpha_est_comp = summary(m1_comp$model)$varcomp['vm(locus, SC, singG = "PSD")', 'component']
Residual_var_comp = summary(m1_comp$model)$varcomp['units!R', 'component']


### Save file ###

# Create a unique stamp for this analysis

unique_stamp = as.character(paste(Sys.info()["nodename"], Sys.time()))
unique_stamp = gsub(" ", "_", unique_stamp)
unique_stamp = gsub(":", "-", unique_stamp)

if(pool_seq){unique_stamp = paste("incorporateQ_pool_seq", pool_seq, "read_length", read_length, "coverage", coverage, "V_logmean", V_logmean, unique_stamp, sep = "_")}

if(verbose){message("Saving data...")}

sim_params = sim_data$sim_params

if(is.null(Ne)){Ne = sim_params$n_ind_exp}
analysis_data = data.frame("proj"=proj, "LDalpha"=LDalpha, "pa"=pa, "Vs"=Vs, "randomise"=randomise, "palpha_method"=palpha, "balpha_method"=paste(balpha[1], balpha[2], sep="_"), "Ne_exp" = paste(Ne, collapse = "_"), "va_true"=parents_info$va_true, "vA_true"=parents_info$vA_true, "vA_est_diag"=vA_est_diag, "vA_est_comp"=vA_est_comp, "palpha_est_diag"=palpha_est_diag, "palpha_est_comp"=palpha_est_comp, "palpha_var_est_diag"=palpha_var_est_diag, "palpha_var_est_comp"=palpha_var_est_comp, "balpha_intercept_est_diag"=balpha_intercept_est_diag, "balpha_intercept_est_comp"=balpha_intercept_est_comp, "balpha_slope_est_diag"=balpha_slope_est_diag, "balpha_slope_est_comp"=balpha_slope_est_comp, "balpha_var_est_diag"=balpha_var_est_diag, "balpha_var_est_comp"=balpha_var_est_comp, "sigma2alpha_est_diag"=sigma2alpha_est_diag, "sigma2alpha_est_comp"=sigma2alpha_est_comp, "Residual_var_diag"=Residual_var_diag, "Residual_var_comp"=Residual_var_comp, "palpha_emp"=parents_info$parameters$palpha, "balpha_intercept_emp"=parents_info$parameters$balpha_0, "balpha_slope_emp"=parents_info$parameters$balpha_1, "sigma2alpha_emp"=parents_info$parameters$sigma2alpha, "analysis_stamp" = unique_stamp)

analysis_data = cbind(sim_params, analysis_data)
write.table(rbind(names(analysis_data), analysis_data), file = paste(output_path, "/Q_diag_complete_", Set_ID, unique_stamp, ".csv", sep = ""),col.names = FALSE, row.names = FALSE, sep = ",")

print(analysis_data)
