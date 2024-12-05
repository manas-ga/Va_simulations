


analyse_parents = function(c_genome,  
                           list_alpha,             # Vector of alphas
                           compute_svdL=FALSE,     # Should SVD of L be computed (overridden if LDalpha=TRUE)
                           LDalpha=FALSE,          # 
                           SNPs,                   # Vector of positions of SNPs
                           RecombRate,             # Recombination rate between adjacent sites
                           HapLength,              # Genome size
                           tol=sqrt(.Machine$double.eps),
                           AtleastOneRecomb=FALSE, # Should there be at least one cross-over per meiosis?
                           calc_nR=TRUE, 
                           verbose=TRUE){ 
  
  ########################
  ## Calculate L and nR ##
  ########################
  
  n0_individuals = nrow(c_genome)/2
  
  if(verbose){
    message("Calculating L...")
  }
  
  paternal<-seq(1, 2*n0_individuals, 2)
  maternal<-paternal+1
  
  c0<-(c_genome[paternal,]+c_genome[maternal,])/2
  
  L<-cov(c0)*(n0_individuals-1)/n0_individuals 
  
  if(compute_svdL | LDalpha){
    if(verbose){
      message("Performing SVD on c0...")
    }  
    svdC<-svd(scale(sqrt(1/n0_individuals)*c0, scale=FALSE), nu = 0)
    retain<-sum(svdC$d>tol)
    UL<-svdC$v[,1:retain]
    DL<-svdC$d[1:retain]
  }
  
  if(calc_nR){
    
    if(verbose){
      message("Calculating the matrix of non-recombinant fractions...")
    }
    
    nR = form_nR(SNPs, RecombRate, HapLength, AtleastOneRecomb)
    
    Lgp<-(cov(c_genome[paternal,])+cov(c_genome[maternal,]))*(n0_individuals-1)/(4*n0_individuals) 
    Ltilde<-Lgp+(1-nR)*(L-Lgp)/nR
    
    rm("Lgp")
    
    
  }else{
    nR = NULL
    Ltilde = NULL
  }
  
  # Segregating sites
  
  seg_sites = length(list_alpha)
  seg_sites_neu = sum(list_alpha==0)
  seg_sites_ben = sum(list_alpha>0)
  seg_sites_del = sum(list_alpha<0)
  
  ##############################################
  ## True additive genic and genetic variance ##
  ##############################################
  
  pbar0 = colMeans(c_genome) 
  diversity = pbar0*(1 - pbar0)/2
  mean_diversity = mean(diversity)
  
  if(verbose){
    message("Calculating the true additive genetic and genic variances...")
  }  
  va_true = sum(diversity*list_alpha^2)          # additive genic variance
  vA_true = t(list_alpha)%*%L%*%list_alpha       # Additive genetic variance
  
  if(verbose){message("Calculating neucleotide diversity (pi)...")}
  
  pair_diff = rep(NA, nrow(c_genome)*(nrow(c_genome) - 1)/2) # vector to record pairwise differences among genomes
  
  count = 1
  
  for(genome1 in 1:(nrow(c_genome)-1)){
    for(genome2 in (genome1+1):nrow(c_genome)){
      pair_diff[count] = sum(abs(c_genome[genome2,] - c_genome[genome1,]))
      count = count + 1 
    }
  }
  
  pi = mean(pair_diff)
  
  if(verbose){message("Computing Watterson's theta...")}
  
  a = sum(1/(1:nrow(c_genome)))
  theta = ncol(c_genome)/a
  
  ###################################################################
  ##### Calculate empirical properties of distribution of alphas ####
  ###################################################################
  
  if(verbose){
    message("Calculating properties of the empirical distribution of alphas...")
  }  
  alpha_properties = alpha_distribution(alpha = list_alpha, p = pbar0)
  
  palpha_emp = alpha_properties$palpha
  balpha_intercept_emp = alpha_properties$balpha_0
  balpha_slope_emp = alpha_properties$balpha_1
  sigma2alpha_emp = alpha_properties$sigma2alpha
  
  # Calculate vA from these empirical properties of alpha
  
  balpha_emp = c(balpha_intercept_emp, balpha_slope_emp)
  
  if(LDalpha){
    TrV<-sum(DL^(2*(palpha_emp+1)))*sigma2alpha_emp
    aLa<-t((alpha_properties$X)%*%balpha_emp)%*%L%*%(alpha_properties$X)%*%balpha_emp-sum(diag(t(alpha_properties$X)%*%L%*%(alpha_properties$X)%*%(alpha_properties$S)))
  }else{
    TrV<-sum(diag(L)^(palpha_emp+1))*sigma2alpha_emp
    aLa<-t((alpha_properties$X)%*%balpha_emp)%*%L%*%(alpha_properties$X)%*%balpha_emp-sum(diag(t(alpha_properties$X)%*%L%*%(alpha_properties$X)%*%(alpha_properties$S)))
  }
  
  vA_alpha_emp<-TrV+aLa
  
  return(list(L=L, Ltilde=Ltilde, nR=nR, svdL=if(compute_svdL | LDalpha){list(UL=UL, DL=DL)}else{NULL}, seg_sites=seg_sites, seg_sites_neu=seg_sites_neu, seg_sites_ben=seg_sites_ben, seg_sites_del=seg_sites_del, mean_diversity=mean_diversity, pi=pi, theta=theta, va_true=va_true, vA_true=vA_true, vA_alpha_emp=vA_alpha_emp, parameters=alpha_properties[1:5]))
  
}