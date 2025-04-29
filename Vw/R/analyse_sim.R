


analyse_sim = function(Set_ID,                # The unique ID of the set of simulations that are controlled by a single R script
                       sim = 1,               # Each set can have multiple sims, but - on the cluster sim must always 1
                       ngen2_optional = NULL, # Allows del_P to be calculated between ngen1 and manually specified ngen2 (which can be different from the last generation)
                       unzip = FALSE,         # Should the SLiM output file be unzipped, read, and then zipped back?
                       slim_output_path,      # The directory where the SLiM outputs (for parents and experimental replicates) are stored (as .txt files)
                       sim_param_path,        # The path to the directory where the .csv file containing simulation parameters is stored
                       extract_genomes_path,  # The path to the python script that extracts genomes and mutations from SLim outputs
                       extract_mut_path,      # The path to the python script that extracts mutations from SLim outputs
                       mutations_path,        # The directory where extracted mutations are to be stored (temp files)
                       c_matrix_path,         # The directory where extracted genomes are to be stored (temp files)
                       output_path,           # The path where the final data file is to be stored
                       n_sample=NULL,         # Number of individuals sampled from the parents' generation (useful if n_ind_exp is large)
                       randomise = TRUE,      # Optionally the reference allele can be randomised
                       delete_temp_files = TRUE,
                       proj = "BLoM",         # projection type for allele frequencies: "LoM", "BLoM", "L" or "N"
                       LDalpha = FALSE,       # Should L or diag(L) be considered while modelling distribution of alphas
                       pa = 1,
                       Vs = "LoNL",           # "L" or "LoNL"
                       method="REML",         # Can be "REML" or "MCMC"
                       palpha = NA,           # If NA palpha is estimated using optim()
                       balpha = c(NA, NA),    # If c(NA,NA) both bedelta intercept and slope are estimated
                       AtleastOneRecomb=FALSE, 
                       Ne=NULL,               # Can be a scalar (same Ne throughout), a vector of length 2 (different Ne's in the neutral (Ne[1]) and selected (Ne[2]) parts of the experiment), or a vector of length ngen2 (different Ne in each generation)
                       predict_Ne = FALSE,    # Should Ne in the experiment be computed using predict_Ne()? If TRUE, provided Ne is replaced by Ne = c(nind_expt, predict_Ne(nind_expt, Ve_w_expt))
                       all.gp = FALSE,        # Ltilde = L'+L''(r/(1-r)) if all.gp=T L'' is assumed 0 and L'=L. 
                       verbose = TRUE
){
  
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
                               delete_temp_files = delete_temp_files)
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
  
  if(predict_Ne){
    Ne = c(sim_data$sim_params$n_ind_exp, predict_Ne(n=sim_data$sim_params$n_ind_exp, Ve=sim_data$sim_params$Ve_w_expt))
  }
  
  ### Fit model ###
  
  if(verbose){message("Performing analyses...")}
  
  m1<-Vw_model(c_genome = NULL,          
               nR = parents_info$nR,
               pbar0 = sim_data$pbar0,                   
               pbar1 = sim_data$pbar1,      
               ngen1=sim_data$ngen1,     
               pbar2 = sim_data$pbar2,       
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
               tol = sqrt(.Machine$double.eps))
  
  
  vA_est = m1$Vw_est 
  palpha_est = m1$palpha 
  palpha_var_est = m1$palpha_var
  balpha_intercept_est = m1$balpha[1]
  balpha_slope_est = m1$balpha[2]
  balpha_var_est = paste(m1$balpha_var[1,1], m1$balpha_var[2,2], m1$balpha_var[1,2], sep = "_")
  sigma2alpha_est = summary(m1$model)$varcomp[1,1]
  
  
  ### Calculate Vw from Buffalo and Coop's method ###
  
  if(sim_data$ngen2-sim_data$ngen1==1){
    
    message("Calculating Vw using Buffalo and Coop's (2019) method (approach 1) ...")
    
    # Three different approaches
    # In each approach either actual (using the approximation of B&C) or exact
    
    ### Approach 1: del_P for all segregating sites and average LD using all segregating sites
    
    BC_fit_1_exact = est_Va_bc(pbar1 = sim_data$pbar1,
                         pbar2 = sim_data$pbar2,
                         L = parents_info$L,
                         nR = parents_info$nR,
                         exact = TRUE)
    
    BC_fit_1 = est_Va_bc(pbar1 = sim_data$pbar1,
                               pbar2 = sim_data$pbar2,
                               L = parents_info$L,
                               nR = parents_info$nR,
                               exact = FALSE)
    
    ### Approach 2: del_P for neutral sites and average LD using all segregating sites
    message("Calculating Vw using Buffalo and Coop's (2019) method (approach 2) ...")
    # Identify selected loci
    
    selected = which(sim_data$list_alpha!=0)
    
    BC_fit_2_exact = est_Va_bc(pbar1 = sim_data$pbar1[,-selected],
                         pbar2 = sim_data$pbar2[,-selected],
                         L = parents_info$L,
                         nR = parents_info$nR,
                         exact = TRUE)
    
    BC_fit_2 = est_Va_bc(pbar1 = sim_data$pbar1[,-selected],
                               pbar2 = sim_data$pbar2[,-selected],
                               L = parents_info$L,
                               nR = parents_info$nR,
                               exact = FALSE)
    
    
    # Approach 3: del_P for neutral sites and average LD using the LD between selected and neutral sites only
    
    message("Calculating Vw using Buffalo and Coop's (2019) method (approach 3) ...")
    BC_fit_3_exact = est_Va_bc(pbar1 = sim_data$pbar1[,-selected],
                         pbar2 = sim_data$pbar2[,-selected],
                         L = parents_info$L,
                         nR = parents_info$nR,
                         selected = selected,
                         exact = TRUE)
    
    BC_fit_3 = est_Va_bc(pbar1 = sim_data$pbar1[,-selected],
                               pbar2 = sim_data$pbar2[,-selected],
                               L = parents_info$L,
                               nR = parents_info$nR,
                               selected = selected,
                               exact = FALSE)
    
    # Combine the results from the three methods separated by "_"
    
    vA_BC = paste("actual", BC_fit_1$vA_BC, BC_fit_2$vA_BC, BC_fit_3$vA_BC,
                  "exact", BC_fit_1_exact$vA_BC, BC_fit_2_exact$vA_BC, BC_fit_3_exact$vA_BC,
                  sep = "_")
    Ne_BC = paste("actual", BC_fit_1$Ne_BC, BC_fit_2$Ne_BC, BC_fit_3$Ne_BC, 
                  "exact", BC_fit_1_exact$Ne_BC, BC_fit_2_exact$Ne_BC, BC_fit_3_exact$Ne_BC, 
                  sep = "_")
  }
  
  ### Save file ###
  
  # Create a unique stamp for this analysis
  
  unique_stamp = as.character(paste(Sys.info()["nodename"], Sys.time()))
  unique_stamp = gsub(" ", "_", unique_stamp)
  unique_stamp = gsub(":", "-", unique_stamp)
  
  if(verbose){message("Saving data...")}
  
  sim_params = sim_data$sim_params
  
  if(is.null(Ne)){Ne = sim_params$n_ind_exp}
  analysis_data = data.frame("proj"=proj, "LDalpha"=LDalpha, "pa"=pa, "Vs"=Vs, "randomise"=randomise, "palpha_method"=palpha, "balpha_method"=paste(balpha[1], balpha[2], sep="_"), "Ne_exp" = paste(Ne, collapse = "_"), "va_true"=parents_info$va_true, "vA_true"=parents_info$vA_true, "vA_est"=vA_est, "vA_alpha_emp"=parents_info$vA_alpha_emp, "vA_BC" = if(sim_data$ngen2-sim_data$ngen1==1){vA_BC}else{NA}, "Ne_BC" = if(sim_data$ngen2-sim_data$ngen1==1){Ne_BC}else{NA}, "Residual_var" = summary(m1$model)$varcomp[2,1], "palpha_emp"=parents_info$parameters$palpha, "balpha_intercept_emp"=parents_info$parameters$balpha_0, "balpha_slope_emp"=parents_info$parameters$balpha_1, "sigma2alpha_emp"=parents_info$parameters$sigma2alpha, "palpha_est"=palpha_est, "palpha_var_est"=palpha_var_est, "balpha_intercept_est"=balpha_intercept_est, "balpha_slope_est"=balpha_slope_est, "balpha_var_est"=balpha_var_est, "sigma2alpha_est"=sigma2alpha_est, "seg_sites"=parents_info$seg_sites, "seg_sites_neu"=parents_info$seg_sites_neu, "seg_sites_ben"=parents_info$seg_sites_ben, "seg_sites_del"=parents_info$seg_sites_del, "mean_diversity"=parents_info$mean_diversity, "theta" = parents_info$theta, "all.gp" = all.gp, "analysis_stamp" = unique_stamp)
  
  analysis_data = cbind(sim_params, analysis_data)
  write.table(rbind(names(analysis_data), analysis_data), file = paste(output_path, "/", Set_ID, "_sim_", sim, "_Data_analysis_", unique_stamp, ".csv", sep = ""),col.names = FALSE, row.names = FALSE, sep = ",")
  
  return(analysis_data)
  
}