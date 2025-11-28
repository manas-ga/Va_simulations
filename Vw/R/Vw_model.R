

Vw_model<-function(c_genome=NULL,    # gamete genotypes (rows gametes (rows 1 & 2 individual 1, rows 3 & 4 individual 2 ...., columns loci) 
                   nR,          # matrix of non-recombinant probabilities between loci
                   pbar0,       # vector of allele frequencies at time-point 0
                   pbar1,       # vector of allele frequencies at time-point 1
                   Q1=NULL,     # list of starting Q matrices for each replicate 
                   ngen1=1,     # number of generations between parents and time-point 1
                   pbar2,       # vector of allele frequencies at time-point 2
                   Q2=NULL,     # list of final Q matrices for each replicate 
                   ngen2,       # number of generations between parents and time-point 2
                   nind=NULL,
                   proj,        # projection type for allele frequencies: "LoM", "BLoM", "L" or "N"
                   LDalpha,
                   pa,
                   palpha,
                   balpha,
                   Vs,
                   method,
                   L,           # between individual covariance in allele proportion
                   Ltilde,      # between gamete covariance in half allele proportion
                   svdL=NULL,   # list with elements UL and DL,
                   NE=NULL,     # Can be a scalar (same Ne throughout), a vector of length 2 (different Ne's in the neutral (Ne[1]) and selected (Ne[2]) parts of the experiment), or a vector of length ngen2 (different Ne in each generation)
                   Ne=NE,
                   diag.projQ=TRUE,
                   tol=sqrt(.Machine$double.eps),
                   save_tprojp=FALSE, 
                   verbose=TRUE)
{
 
  if(is.null(nind)){
    if(is.null(c_genome)){
      stop("nind must be specified if c_genome is null")
    }else{
      nind<-nrow(c_genome)/2
    }
  }else{
    if(!is.null(c_genome)){
      if(nind!=nrow(c_genome)/2){
        stop("nind is not equal to half the number of rows in c_genome ")
      }
    }  
  } 

  if(is.null(NE)){NE = nind}   # If Ne is not provided, Ne should be nind

  if(is.null(Ne)){Ne = NE}
  
  if(length(NE)!=length(Ne)){stop("NE and Ne should be the same length")}

  if(!length(Ne)%in%c(1, 2, ngen2)){stop("Ne must be a vector of length either 1, 2, or ngen2")}
  
  asreml.options(Cfixed = TRUE)
  
  if(is.null(L) | is.null(Ltilde) | (is.null(svdL) & (proj=="BLoM" | LDalpha))){
    
    if(is.null(c_genome)){
      stop("c_genome is required if L, Ltilde or svdL are NULL")
    }
    
    if(is.null(L) | (is.null(svdL) & (proj=="BLoM" | LDalpha))){
      
      paternal<-seq(1, 2*nind, 2)
      maternal<-paternal+1
      
      c0<-(c_genome[paternal,]+c_genome[maternal,])/2
      
      if(is.null(L)){ 
        if(verbose){
          message("Computing L in the parents' generation...")
        }
        L<-cov(c0)*(nind-1)/nind
      }  
    }
    
    if(is.null(Ltilde)){
      Lgp<-(cov(c_genome[paternal,])+cov(c_genome[maternal,]))*(nind-1)/(4*nind) 
      Ltilde<-Lgp+(1-nR)*(L-Lgp)/nR
      rm("Lgp")
    }
  }
  
  
  if(!proj%in%c("LoM", "BLoM", "L", "N")){stop("proj must be one of 'LoM', 'L', 'N'")}
  if(!Vs%in%c("LoNL", "L")){stop("Vs must be either 'LoNL' or 'L'")}
  if(!method%in%c("REML", "MCMC")){stop("method must be either 'REML' or 'MCMC'")}
  
  #####################################
  ### calculate projection matrices ###
  #####################################
  
  nsnps<-ncol(pbar1)
  nrep<-nrow(pbar1)
  
  if(ncol(pbar1)!=nsnps){stop("pbar1 should have as many columns as the number of columns in Co (the number of SNPs)")}
  if(ncol(pbar2)!=nsnps){stop("pbar2 should have as many columns as the number of columns in Co (the number of SNPs)")}
  
  if(nrow(pbar1)!=nrow(pbar2)){stop("pbar1 and pbar2 should have the same number of rows (the number of replicates)")}
  
  if(proj=="BLoM" | LDalpha){  # singular vectors of L required 
    
    if(is.null(svdL)){
      
      if(verbose){
        message("Performing SVD on C0...")
      }
      
      c0=scale(sqrt(1/nind)*c0, scale=FALSE)
      
      svdC<-svd(c0, nu = 0)
      retain<-sum(svdC$d>tol)
      UL<-svdC$v[,1:retain]
      DL<-svdC$d[1:retain]
      rm("svdC")
      gc(verbose = FALSE)
      
    }else{
      UL<-svdL$UL
      DL<-svdL$DL
      rm("svdL")
      gc(verbose = FALSE)
      retain<-ncol(UL)
    }
  }
  
  ##################################
  ### Calculate Drift covariance ###
  ##################################
  
  if(proj=="LoM" | proj=="BLoM"){ 
    
    if(verbose){
      message("Computing M...")
    }
    
    # Calculate the summation in two steps to save memory: (1) Write down the first term, (2) add the remaining terms only if ngen2 - ngen1 >1  Note this works even if ngen1=0, at least when Ne is constant.
    
    ### If Ne is a scalar, the same Ne is to be used throughout ###
    
    if (length(Ne) == 1){
            
      M = ((1-1/(2*Ne))^ngen1)*(1/NE)*nR^(1+ngen1)
      
      # Only perform further summations if (ngen2-ngen1 > 1)
      
      if((ngen2 - ngen1) > 1){
        for (x in (ngen1 + 1):(ngen2 - 1)){
          M = M + ((1-1/(2*Ne))^x)*(1/NE)*nR^(x+1)
        }       
      }
      
    }
    
    if(length(Ne)==2){
      
      M = ((1-1/(2*Ne[1]))^ngen1)*(1/NE[2])*nR^(1+ngen1)
      
      # Only perform further summations if (ngen2-ngen1 > 1)
      
      if((ngen2 - ngen1) > 1){
        for (x in (ngen1 + 1):(ngen2 - 1)){
          M = M + ((1-1/(2*Ne[1]))^ngen1)*((1-1/(2*Ne[2]))^(x-ngen1))*(1/NE[2])*nR^(x+1)
        }       
      }
      
      
    }
    
    if(length(Ne)==ngen2){
      
      
      M = (prod(1-1/(2*Ne[1:ngen1])))*(1/NE[ngen1+1])*nR^(1+ngen1)
           
      # Only perform further summations if (ngen2-ngen1 > 1)
     
      if((ngen2 - ngen1) > 1){
        for (x in (ngen1 + 1):(ngen2 - 1)){
          M = M + (prod(1-1/(2*Ne[1:x])))*(1/NE[x+1])*nR^(1+x)
        }       
     }     
      
    }
 
    Drift<-Ltilde*M
    rm("M")
    # Garbage collection
    gc(verbose = FALSE)
  }
  
  ######################################
  ### Calculate Selec (\mathcal{L})  ###
  ######################################
  
  if(Vs=="LoNL"){ 
    
    if(verbose){
      message("Computing N...")
    }
    
    # Calculate the summation in two steps to save memory: (1) Write down the first term, (2) add the remaining terms only if ngen2 - ngen1 >1
    
    if(length(Ne)==1){
      
      
      if(ngen1!=0){  
        N = {((1-1/(2*Ne))^ngen1)*nR^ngen1}
      }else{
        if(ngen2>1){
          N = matrix(0, ncol(nR), ncol(nR))
        }
      }
      
      # Only perform further summations if (ngen2-ngen1 > 1)
      
      if((ngen2 - ngen1) > 1){
        
        for (x in (ngen1 + 1):(ngen2-1)){
          N = N + ((1-1/(2*Ne))^x)*nR^x
        }
        
      }
      
    }
    
    if(length(Ne)==2){
      
      if(ngen1!=0){  
        N = {((1-1/(2*Ne[1]))^ngen1)*nR^(ngen1)}
      }else{
        if(ngen2>1){
          N = matrix(0, ncol(nR), ncol(nR))
        }
      }
      
      # Only perform further summations if (ngen2-ngen1 > 1)
      
      if((ngen2 - ngen1) > 1){
        for (x in (ngen1 + 1):(ngen2-1)){
          N = N + ((1-1/(2*Ne[1]))^ngen1)*((1-1/(2*Ne[2]))^(x-ngen1))*nR^x
          
        }
        
      }
      
    }
    
    if(length(Ne)==ngen2){
      
      if(ngen1!=0){  
        N = (prod(1-1/(2*Ne[1:ngen1])))*nR^ngen1
      }else{
        if(ngen2>1){
          N = matrix(0, ncol(nR), ncol(nR))
        }
      }
      
      # Only perform further summations if (ngen2-ngen1 > 1)
      
      if((ngen2 - ngen1) > 1){
        
        for (x in (ngen1 + 1):(ngen2-1)){
          N = N + (prod(1-1/(2*Ne[1:x])))*nR^x
        }
      }
      
    }
    
    if(ngen1==0){
      Selec<-L
      if(ngen2>1){
        Selec<-Selec+Ltilde*N
        rm("N")
        gc(verbose = FALSE)   
      }
    }else{
      Selec<-Ltilde*N
      rm("N")
      gc(verbose = FALSE)   
    }
    
  }else{
    Selec<-L*(ngen2-ngen1)
  }  
  
  if(proj=="LoM"){
    
    if(verbose){
      message("Performing eigendecomposition of L*M...")
    }
    
    sdLoM<-RSpectra::eigs(Drift, min(nind, nsnps))
    retain<-sum(sqrt(sdLoM$values)>tol)
    U<-sdLoM$vectors[,1:retain]
    D<-sqrt(sdLoM$values[1:retain])
    rm("sdLoM")
    # Garbage collection
    gc(verbose = FALSE)
  }
  
  if(proj=="BLoM"){ 
    
    if(verbose){
      message("Performing eigendecomposition of t(U)%*%(L*M)%*%U...")
    }
    
    sdLoM<-eigen(t(UL)%*%Drift%*%UL)
    retain<-sum(sqrt(sdLoM$values)>tol)
    U<-UL%*%sdLoM$vectors[,1:retain]
    D<-sqrt(sdLoM$values[1:retain])
    rm("sdLoM")
    # Garbage collection
    gc(verbose = FALSE)
    
  }   
  
  if(verbose){
    message("Computing the projection matrix...")
  }
  
  if(proj=="N"){
    
    tprojp<-diag(nrow(pbar1))
    
  }else{
    
    tprojp<-U%*%diag(D^(-pa))
    
  }
  
  rm(list = c("U", "D"))
  # Garbage collection
  gc(verbose = FALSE)
  
  if(!is.null(Q1)){

    for(i in 1:length(Q1)){
       if(class(Q1[[i]])=="ddiMatrix"){
         Q1[[i]]<-t(tprojp)%*%Diagonal(nrow(Q1[[i]]), 2*diag(Ltilde)*(Q1[[i]]@x+Q2[[i]]@x))%*%tprojp 
       }else{
         Q1[[i]]<-2*t(tprojp)%*%(Ltilde*(Q1[[i]]+Q2[[i]]))%*%tprojp
       }
    }
    projQ<-bdiag(Q1)
    rm(list = c("Q1", "Q2"))
    gc(verbose = FALSE)
    
    if(diag.projQ){
      projQ<-Diagonal(nrow(projQ), diag(projQ))
    }
  }else{
    projQ<-NULL
  }
  
  if(is.na(palpha)){
    
    if(verbose){
      message("Estimating palpha...")
    }
    
    palpha<-optim(0, fit.model, balpha=balpha, LDalpha = LDalpha, nsnps=nsnps, UL=UL, DL=DL, L=L, ngen2=ngen2, ngen1=ngen1, tprojp=tprojp, pbar0=pbar0, pbar1=pbar1, pbar2=pbar2, projQ=projQ, nrep=nrep, LLonly=TRUE, Selec=Selec, verbose=verbose, method = "L-BFGS-B", lower = -2, upper =2, control = list(fnscale=-1, factr = 1e+11), hessian=TRUE)
    
    palpha_var<--1/palpha$hessian
    palpha<-palpha$par
    
  }else{
    palpha_var<-0
  }
  
  if(verbose){
    message("Fitting the final model...")
  }
  
  output<-fit.model(palpha=palpha, balpha=balpha, LDalpha = LDalpha, nsnps=nsnps, UL=UL, DL=DL, L=L, ngen2=ngen2, ngen1=ngen1, tprojp=tprojp, pbar0=pbar0, pbar1=pbar1, pbar2=pbar2, projQ=projQ, nrep=nrep, LLonly=FALSE, Selec=Selec, verbose=verbose)
  
  if(verbose){
    message("Calculating the estimate of Vw...")
  }
  
  if(method=="REML"){
    
    sigma2alpha<-summary(output$model)$varcomp['vm(locus, SC, singG = "PSD")',1]/output$SC_scale

    S<-matrix(0,2,2)
    colnames(S)<-rownames(S)<-c("int", "pmq")

    if(is.na(balpha[1]) & is.na(balpha[2])){
      balpha<-summary(output$model, coef=TRUE)$coef.fixed[c("int", "pmq"),1]
      S<-output$model$Cfixed[c("int", "pmq"),c("int", "pmq")]
    }
    if(is.na(balpha[1]) & !is.na(balpha[2])){
      balpha[1]<-summary(output$model, coef=TRUE)$coef.fixed["int",1]
      S[1,1]<-output$model$Cfixed["int", "int"]
    }
    if(!is.na(balpha[1]) & is.na(balpha[2])){
      balpha[2]<-summary(output$model, coef=TRUE)$coef.fixed["pmq",1]
      S[2,2]<-output$model$Cfixed["pmq", "pmq"]
    }
  }
  if(method=="MCMC"){
    sigma2alpha<-posterior.mode(output$model$VCV[,1])
    balpha<-colMeans(output$model$Sol[,c("int", "pmq")])
  }
  
  X<-cbind(rep(1, length(pbar0)), 2*pbar0-1)
  
  if(LDalpha){
    TrV<-sum(DL^(2*(output$palpha+1)))*sigma2alpha
    aLa<-t(X%*%balpha)%*%L%*%X%*%balpha-sum(diag((t(X)%*%L)%*%(X%*%S)))
  }else{
    TrV<-sum(diag(L)^(output$palpha+1))*sigma2alpha
    aLa<-t(X%*%balpha)%*%L%*%X%*%balpha-sum(diag((t(X)%*%L)%*%(X%*%S)))
    
  }
  
  Vw_est<-TrV+aLa
  
  
  
  return(list(Vw_est=Vw_est, data=output$data, model=output$model, SC=output$SC, palpha=output$palpha, balpha=balpha, palpha_var=palpha_var, balpha_var=S, tprojp=if(save_tprojp){tprojp}else{NULL}, X=X, DL=ifelse(exists("DL"), DL, NA)))
  
}