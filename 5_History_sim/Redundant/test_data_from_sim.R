### Base path and path to Vw.Rmd (file containing Jarrod's functions) (depending on the system) ###

if(Sys.info()["nodename"]%in%c("bigfoot", "bigshot", "bigbird", "bigyin", "biggar", "bigwig", "c1", "c2", "c3", "c4", "c5", "c6")){
  
  base_path = "/ceph/users/marun/Va_simulations/5_History_sim"
  Vw_path = "/ceph/users/marun/Va_simulations/6_Code_Test/Vw.Rmd"
  
}else{
  
  
  if(Sys.info()["nodename"]=="vera.bio.ed.ac.uk"){
    
    base_path = "/data/home/msamant/Manas/Va_simulations/Github/Va_simulations/5_History_sim" ## ON VERA
    Vw_path = "/data/home/msamant/Manas/Va_simulations/Github/Va_simulations/6_Code_Test/Vw.Rmd"  ### Jarrod's functions and other code is stored here
    
    
  }else{
    
    if(Sys.info()["sysname"]=="Linux"){
      
      base_path = "/mnt/c/Users/msamant/Documents/GitHub/Va_simulations/5_History_sim" ## Local Wsl
      Vw_path = "/mnt/c/Users/msamant/Documents/GitHub/Va_simulations/6_Code_test/Vw.Rmd" ### Jarrod's functions and other code is stored here
      
    }else{
      
      base_path = "C:/Users/msamant/Documents/GitHub/Va_simulations/5_History_sim" ## Local windows
      Vw_path = "C:/Users/msamant/Documents/GitHub/Va_simulations/6_Code_test/Vw.Rmd" ### Jarrod's functions and other code is stored here
      
    }
    
  }
  
  
}

####################################
######### Packages #################
####################################

library(MCMCglmm)
library(asreml)
library(Matrix)
library(rmutil)
library(pryr) ## For tracking memory usage using mem_used()
library(bigalgebra)

#################################
#### Load Jarrod's functions ####
#################################

functions_only=TRUE ## Read only the functions

rmarkdown::render(file.path(Vw_path))

# Load the mutations in the parents' generation

mutations_0 = read.table(file.choose(), header = T)
mutations_0 = mutations_0[order(mutations_0$Temp_ID),]
list_alpha = 2*(mutations_0$s)
p_parents = mutations_0$Number/2000 # Allele frequencies in parents

# Randomise reference alleles 

ran_vect = sample(c(0, -1), nrow(mutations_0),  replace = T)
ran_vect_L = ifelse(ran_vect==-1, -1, 1)
list_alpha = list_alpha*ran_vect_L
p_parents = abs(p_parents + ran_vect)

alpha_properties = alpha_distribution(alpha = list_alpha, p = p_parents)

pdelta_emp = alpha_properties$pdelta
bdelta_intercept_emp = alpha_properties$bdelta_int
bdelta_slope_emp = alpha_properties$bdelta_slope
sigma2delta_emp = alpha_properties$sigma2delta

# Load the final data

d = read.csv(file.choose(), header=T)
d = d[grep("Set_ID", d$Set_ID, invert=TRUE),]

d$vA_est = as.numeric(d$vA_est)
d$vA_true = as.numeric(d$vA_true)
d$vA_left = as.numeric(d$vA_true)
d$pdelta_est = as.numeric(d$pdelta_est)
d$bdelta_intercept_est = as.numeric(d$bdelta_intercept_est)
d$bdelta_slope_est = as.numeric(d$bdelta_slope_est)
d$sigma2delta_est = as.numeric(d$sigma2delta_est)
d$r_expt = as.numeric(d$r_expt)
d$r = as.numeric(d$r)
d$sequence_length = as.numeric(d$sequence_length)
d$mu = as.numeric(d$mu)
d$s_pmq = as.numeric(d$s_pmq)
d$seg_sites = as.numeric(d$seg_sites)
d$seg_sites_neu = as.numeric(d$seg_sites_neu)
d$seg_sites_del = as.numeric(d$seg_sites_del)
d$seg_sites_ben = as.numeric(d$seg_sites_ben)
d$pdelta_emp = as.numeric(d$pdelta_emp)
d$bdelta_intercept_emp = as.numeric(d$bdelta_intercept_emp)
d$bdelta_slope_emp = as.numeric(d$bdelta_slope_emp)
d$sigma2delta_emp = as.numeric(d$sigma2delta_emp)
d$s_pmq = as.numeric(d$s_pmq)
d$n_cages = factor(d$n_cages, levels = c("3", "5", "10"))
d$n_ind_exp = factor(d$n_ind_exp, levels = c("100", "500", "1000"))
d$ngen_expt = factor(d$ngen_expt, levels = c("1", "3", "5"))
d$map_length_history = factor(d$r*d$sequence_length)
d$map_length_expt = factor(d$r_expt*d$sequence_length)

parameters_est = list("pdelta"=d$pdelta_est, "bdelta_int" = d$bdelta_intercept_est, "bdelta_slope" = d$bdelta_slope_est, "sigma2delta" = d$sigma2delta_est)

plot_alpha_distribution(alpha = list_alpha, p = p_parents, parameters = alpha_properties,pch=16, cex=0.1, col="grey")
plot_alpha_distribution(alpha = list_alpha, p = p_parents, parameters = parameters_est, pch=16, cex=0.2, col="grey")
