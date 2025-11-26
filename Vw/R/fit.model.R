
fit.model<-function(palpha, balpha, LDalpha, nsnps, UL, DL, L, ngen2, ngen1, tprojp, pbar0, pbar1, pbar2, projQ=NULL, nrep, Selec, LLonly=FALSE, method = "REML", verbose=TRUE){
  
  if(verbose){
    message("Computing the covariance structure of locus effects...")
    message("Computing covp over 1 generation...")
  }
  
  if(LDalpha){
    if(palpha==0){
      covp<-diag(nsnps)
    }else{
      covp<-UL%*%diag(DL^(2*palpha))%*%t(UL)
    }  
  }else{
    covp<-diag(L)^palpha
  }
  
  rm(list = c("UL", "DL"))
  gc(verbose = FALSE)
  pmq<-2*pbar0-1
  
  pmq<-Selec%*%pmq
  int<-Selec%*%rep(1, nsnps)
  
  if(verbose){
    message("Computing SC...")
  }
  
  projpSelec<-t(tprojp)%*%Selec
  
  rm("Selec")
  gc(verbose = FALSE)
  
  if(LDalpha){
    SC<-projpSelec%*%covp
    
    rm("covp")
    gc(verbose = FALSE)
    
    SC<-SC%*%t(projpSelec)
    # computing COV(\Delta p_m, \Delta p_n) on the projected scale.
  }else{
    SC<-projpSelec%*%(t(projpSelec)*covp)
    # if LDalpha=F, covp is diagonal and this is more efficient. 
  }
  SC_scale<-mean(diag(SC))
  SC<-SC/SC_scale

  attr(SC, "INVERSE")<-FALSE
  dimnames(SC) <- list(1:nrow(SC),1:nrow(SC))  # used for full-form matrices

  if(!is.null(projQ)){
    projQ<-as(projQ, "TsparseMatrix")
    attr(projQ, "INVERSE")<-FALSE
    dimnames(projQ) <- list(1:nrow(projQ),1:nrow(projQ))
  } 

  if(verbose){
    message("Projecting allele frequencies...")
  }
  
  # Before projecting allele frequencies, the positions where either pbar1 or pbar2 has an NA, need to be turned into 0s in both pbar1 and pbar2
  # Identify the positions to be turned into zeros
  
  zeros = is.na(pbar1*pbar2)
  pbar1[zeros] = 0
  pbar2[zeros] = 0
  
  pbar1_proj<-pbar1%*%tprojp
  pbar2_proj<-pbar2%*%tprojp
  pmq_proj<-t(matrix(pmq, nsnps,nrep))%*%tprojp
  int_proj<-t(matrix(int, nsnps,nrep))%*%tprojp
  
 # dat.gaussian<-data.frame(delta=c(pbar2_proj-pbar1_proj), locus=gl(ncol(pbar1_proj),nrep,ncol(pbar1_proj)*nrep), rep=gl(nrep,1,ncol(pbar1_proj)*nrep), pmq=c(pmq_proj), int=c(int_proj))
  
  dat.gaussian<-data.frame(delta=c(t(pbar2_proj-pbar1_proj)), locus=gl(ncol(pbar1_proj),1,ncol(pbar1_proj)*nrep), rep=gl(nrep,ncol(pbar1_proj),ncol(pbar1_proj)*nrep), pmq=c(t(pmq_proj)), int=c(t(int_proj)))

  
  if(!is.na(balpha[1])){
    dat.gaussian$int<-dat.gaussian$int*balpha[1]
  }
  if(!is.na(balpha[2])){
    dat.gaussian$pmq<-dat.gaussian$pmq*balpha[2]
  }
  
  ##############
  # Fit models #
  ##############
  
  if(verbose){
    message("Fitting linear mixed model...")
  }
  
  prior<-list(R=list(V=1, nu=0), G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)))
  
  SC[lower.tri(SC)]<-t(SC)[lower.tri(SC)]
  # force symmetry in case it is not exactly symmetric
  
  if(method=="REML"){

    if(is.null(projQ)){
      random = ~vm(locus, SC, singG="PSD")
    }else{
      random = ~vm(locus, SC, singG="PSD")+vm(units, projQ, singG="PSD")
    }
 
    if(is.na(balpha[1]) & is.na(balpha[2])){
      m1<-asreml(delta~int+pmq-1, random = random, data=dat.gaussian, Cfixed=TRUE)
    }
    if(is.na(balpha[1]) & !is.na(balpha[2])){
      m1<-asreml(delta~int+offset(pmq)-1, random = random, data=dat.gaussian, Cfixed=TRUE)
    }
    if(!is.na(balpha[1]) & is.na(balpha[2])){
      if(balpha[1]==0){
        m1<-asreml(delta~pmq-1, random = random, data=dat.gaussian, Cfixed=TRUE)
      }else{
        m1<-asreml(delta~offset(int)+pmq-1, random = random, data=dat.gaussian, Cfixed=TRUE)
      }
    }  
    if(!is.na(balpha[1]) & !is.na(balpha[2])){
      warning("asreml doesn't allow models without fixed effects, so intercept fitted but replaced with balpha[1]!")
      m1<-asreml(delta~offset(pmq + int), random = random, data=dat.gaussian, Cfixed=TRUE)
    } 
  }
  
  if(method=="MCMC"){
    
    if(LLonly){stop("method = MCMC specified so can't return log-likelihood with LLony=TRUE")}
    
    invSC<-solve(SC)
    invSC <- as(invSC, "sparseMatrix") 
    attr(invSC, "rowNames") <- 1:ncol(SC)
    attr(invSC, "colNames") <- 1:ncol(SC)
    
    prior<-list(B=list(mu=balpha, V=diag(2)*1e+10))
    diag(prior$B$V)[which(!is.na(balpha))]<-1e-10
    prior$B$mu[which(is.na(balpha))]<-0
    
    m1<-MCMCglmm(delta~pmq+int-1, random=~locus, data=dat.gaussian, ginverse=list(locus=invSC), family="gaussian", pr=TRUE, prior=prior)
  }
  
  if(LLonly){
    return(m1$loglik)
  }else{
    return(list(data=dat.gaussian, model=m1, SC=SC, palpha=palpha, SC_scale=SC_scale))
  }
}
