
sim_exp<-function(c_genome, nind, r, position, nrep, alpha=NULL, Ve=0, fitness_model="Exp"){

  if(ncol(c_genome)!=length(position)){stop("position should be the same length as number of columns in c_genome")}

  if(!is.null(alpha)){
    if(ncol(c_genome)!=length(alpha)){stop("alpha should be the same length as number of columns in c_genome")}
  }

  reordering<-order(position)

  c_genome<-c_genome[,reordering]
  position<-position[reordering]

  # reorder SNPs according to map position

  n0_individuals<-nrow(c_genome)/2
  n_loci<-ncol(c_genome)

  pbar<-matrix(NA, nrep, n_loci)

  if(!is.null(alpha)){
    alpha<-alpha[reordering]
    genetic_fitness<-c_genome%*%alpha/2
    genetic_fitness<-genetic_fitness[seq(1, 2*n0_individuals, 2)]+genetic_fitness[seq(2, 2*n0_individuals, 2)]
  }else{
    genetic_fitness<-rep(0, n0_individuals)
  }

  genetic_fitness<-1+genetic_fitness-mean(genetic_fitness)
  # rescale to have mean of 1 (not necessary under the exp model)  

  if(any(genetic_fitness<0)){
    warning(paste(sum(genetic_fitness<0), "genetic fitnesses negative - set to zero"))
    genetic_fitness[which(genetic_fitness<0)]<-0
  }

  for(rep in 1:nrep){

    fitness<-genetic_fitness
    if(Ve!=0){
      if(fitness_model=="Exp"){
          fitness<-fitness+rnorm(n0_individuals, 0, sqrt(Ve))
      }
      if(fitness_model=="I"){
          m<-min(fitness)
          fitness<-fitness+rgamma(n0_individuals, shape=(m^2)/Ve, scale=Ve/m)-m
          # We wish to generate random environmental variables, X, with E[X]=0,  VAR(X)=Ve and min(X) = -min(genetic_fitness) =-m  so that fitness is non-negative 

          # To do so we simulate Y from a gamma and generate X = Y-m to satisfy the constraint that X+G>0 where G is fitness based on genotype.

          # Have a=shape and s=scale then E[X] = a*s-m and VAR(X) = a*s^2

          # Since E[X]=0 then a*s = m and so VAR(X) = m*s and s=VAR(X)/m=Ve/m

          # Since a*s = m then a*Ve/m = m and so a = (m^2)/Ve
      }
    }

    if(fitness_model=="Exp"){
      fitness<-exp(fitness)
    }
      
    genomes<-matrix(0, nind, n_loci)

    for(i in 1:nind){

      for(j in 1:2){  # sample gametes from both parents

        parent<-sample(1:n0_individuals, 1, prob=fitness)
        # sample parent j of individual i

        gamete<-rep(rbinom(1, prob=0.5, size=1), n_loci)
        # sample grandparental gamete of parent j of individual i

        r_event<-reda::simEventData(rho=r, end=max(position))
        # simulate recombination events

        if(nrow(r_event)>0){
          r_event<-r_event$time[which(r_event$event==1)]
          for(k in 1:length(r_event)){
            gamete[which(position>r_event[k])]<-abs(1-gamete[which(position>r_event[k])])
          } 
        }
        # switch grandparental gamete when a recombination event happens and so grandparental gamete is a string of 0's and 1's 

        gam0<-which(gamete==0)
        # which positions does individual i inherit from grandparental gamete 0 of parent j

        if(length(gam0)!=0){     
          genomes[i,gam0]<-genomes[i,gam0]+c_genome[(parent-1)*2+1,gam0]
          # add any reference alleles inherited from grandparental gamete 0 
        }
        if(length(gam0)!=n_loci){
          if(length(gam0)!=0){
            genomes[i,-gam0]<-genomes[i,-gam0]+c_genome[(parent-1)*2+2,-gam0]
          }else{
            genomes[i,]<-genomes[i,]+c_genome[(parent-1)*2+2,]
          }  
          # add any reference alleles inherited from grandparental gamete 1
        }
      }
    }  
    pbar[rep,]<-colMeans(genomes)/2
  }  
  pbar<-pbar[,order(reordering)]
  return(pbar)
}


predict_Ne<-function(Ve, n, fitness_model = "Exp"){

  if(fitness_model == "Exp"){
    # assuming meanlog=0
    Vo<-4*exp(Ve)-2
  }
  if(fitness_model == "I"){
    # assuming mean=1
    Vo<-4*Ve+2
  }
  Ne<-4*n/(2+Vo)
  return(Ne)
}

form_nR<-function(SNPs,                   # SNP positions
                  RecombRate,             # Recombination rate between adjacent sites
                  HapLength,              # Genome size
                  AtleastOneRecomb        # Should there be at least one cross-over per meiosis?
                  ){

  nsnps<-length(SNPs)
  
  message("Computing the distance matrix...")

  Dist<-as.matrix(dist(SNPs, diag=TRUE, upper=TRUE))
  # conditional on the total number of recombination events (tn) along a chromosome, the number that fall in the interval is binomial with probability equal to the relative length of the interval to the whole and number of trials equal to tn. The non-recombination fraction is then the probability that an even number of events (including zero) fall in the interval. 

  message("Computing nR...")
  
  if(AtleastOneRecomb){

    maxn<-actuar::qztpois(1-1e-5, RecombRate*HapLength)
    # maximum number of total recombination events likely to be seen
  
    nR<-matrix(0, nsnps, nsnps)

    for(nt in 0:maxn){
      # could iterate nt from 1 since Pr(0 recombination events in total)=0.
      for(ne in seq(0,maxn,2)){
         nR<-nR+dbinom(ne, prob=Dist/HapLength, size=nt)*actuar::dztpois(nt, RecombRate*HapLength)
         # Pr(even number of crossovers falling between snps given nt)*Pr(nt recombination events in total) then sum over nt
      }
    }
    nR<-nR/actuar::pztpois(maxn, RecombRate*HapLength)
    # rescale so diag(R)=1
  
  }else{

    nR<-1 - 0.5*(1 - exp(-2*Dist*RecombRate))
  }

  return(nR)

}

fit.model<-function(palpha, balpha, LDalpha, nsnps, UL, DL, L, ngen2, ngen1, nind, tprojp, pbar0, pbar1, pbar2, nrep, Selec, Ne_factor=1, LLonly=FALSE, method = "REML", verbose=TRUE){
  
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

  attr(SC, "INVERSE")<-FALSE
  dimnames(SC) <- list(1:nrow(SC),1:nrow(SC))  # used for full-form matrices
  
  if(verbose){
    message("Projecting allele frequencies...")
  }

  pbar1_proj<-pbar1%*%tprojp
  pbar2_proj<-pbar2%*%tprojp
  pmq_proj<-t(matrix(pmq, nsnps,nrep))%*%tprojp
  int_proj<-t(matrix(int, nsnps,nrep))%*%tprojp

  dat.gaussian<-data.frame(delta=c(pbar2_proj-pbar1_proj), locus=gl(ncol(pbar1_proj),nrep,ncol(pbar1_proj)*nrep), rep=gl(nrep,1,ncol(pbar1_proj)*nrep), pmq=c(pmq_proj), int=c(int_proj))


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
    if(is.na(balpha[1]) & is.na(balpha[2])){
        m1<-asreml(delta~int+pmq-1, random = ~vm(locus, SC, singG="PSD"), data=dat.gaussian, Cfixed=TRUE)
    }
    if(is.na(balpha[1]) & !is.na(balpha[2])){
        m1<-asreml(delta~int+offset(pmq)-1, random = ~vm(locus, SC, singG="PSD"), data=dat.gaussian, Cfixed=TRUE)
    }
    if(!is.na(balpha[1]) & is.na(balpha[2])){
        m1<-asreml(delta~offset(int)+pmq-1, random = ~vm(locus, SC, singG="PSD"), data=dat.gaussian, Cfixed=TRUE)
    }
    if(!is.na(balpha[1]) & !is.na(balpha[2])){
      warning("asreml doesn't allow models without fixed effects, so intercept fitted but replaced with balpha[1]!")
      m1<-asreml(delta~offset(pmq + int), random = ~vm(locus, SC, singG="PSD"), data=dat.gaussian, Cfixed=TRUE)
    } 
  }

  if(method=="MCMC"){

    if(LLonly){stop("method = MCMC specified so can't return log-likelihood with LLony=TRUE")}

    invSC<-solve(SC)
    invSC <- as(invSC, "sparseMatrix") 
    attr(invSC, "rowNames") <- 1:retain
    attr(invSC, "colNames") <- 1:retain

    prior<-list(B=list(mu=balpha, V=diag(2)*1e+10))
    diag(prior$B$V)[which(!is.na(balpha))]<-1e-10
    prior$B$mu[which(is.na(balpha))]<-0

    m1<-MCMCglmm(delta~pmq+int-1, random=~locus, data=dat.gaussian, ginverse=list(locus=invSC), family="gaussian", pr=TRUE, prior=prior)
  }
  
  if(LLonly){
    return(m1$loglik)
  }else{
    return(list(data=dat.gaussian, model=m1, SC=SC, palpha=palpha))
  }
}

Vw_model<-function(c_genome,    # gamete genotypes (rows gametes (rows 1 & 2 individual 1, rows 3 & 4 individual 2 ...., columns loci) 
                   nR,          # matrix of non-recombinant probabilities between loci
                   pbar0,       # vector of allele frequencies at time-point 0
                   pbar1,       # vector of allele frequencies at time-point 1
                   ngen1=1,     # number of generations between parents and time-point 1
                   pbar2,       # vector of allele frequencies at time-point 2
                   ngen2,       # number of generations between parents and time-point 2
                   nind,        # population size in each replicate
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
                   Ne_factor=1,
                   tol=sqrt(.Machine$double.eps),
                   save_tprojp=FALSE, 
                   verbose=TRUE)
  {
  
  asreml.options(Cfixed = TRUE)

  if(is.null(L) | is.null(Ltilde) | (is.null(svdL) & (proj=="BLoM" | LDalpha))){

    if(is.null(c_genome)){
      stop("c_genome is required if L, Ltilde or svdL are NULL")
    }

    n0_individuals<-nrow(c_genome)/2

    if(is.null(L) | (is.null(svdL) & (proj=="BLoM" | LDalpha))){

      paternal<-seq(1, 2*n0_individuals, 2)
      maternal<-paternal+1

      c0<-(c_genome[paternal,]+c_genome[maternal,])/2

      if(is.null(L)){ 
        if(verbose){
          message("Computing L in the parents' generation...")
        }
        L<-cov(c0)*(n0_individuals-1)/n0_individuals
      }  
    }

    if(is.null(Ltilde)){
      Lgp<-(cov(c_genome[paternal,])+cov(c_genome[maternal,]))*(n0_individuals-1)/(4*n0_individuals) 
      Ltilde<-Lgp+(1-nR)*(L-Lgp)/nR
      rm("Lgp")
    }
  }



  if(!proj%in%c("LoM", "BLoM", "L", "N")){stop("proj must be one of 'LoM', 'L', 'N'")}
  if(!Vs%in%c("LoNL", "L")){stop("Vs must be either 'LoNL' or 'L'")}
  if(!method%in%c("REML", "MCMC")){stop("method must be either 'REML' or 'MCMC'")}

  #################################
  # calculate projection matrices #
  #################################
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

  if(proj=="LoM" | proj=="BLoM"){ 
    
     if(verbose){
       message("Computing M...")
     }

     # Calculate the summation in two steps to save memory: (1) Write down the first term, (2) add the remaining terms only if ngen2 - ngen1 >1  Note this works even if ngen1=0, at least when Ne is constant.

     M = ((1-1/(2*nind*Ne_factor))^(ngen1))*(1/(nind*Ne_factor))*nR^(1+ngen1)
     
     # Only perform further summations if (ngen2-ngen1 > 1)
     
     if((ngen2 - ngen1) > 1){
       for (x in (ngen1 + 2):(ngen2)){
        M = M + ((1-1/(2*nind*Ne_factor))^(x-1))*(1/(nind*Ne_factor))*nR^x
      }       
     }

     Drift<-Ltilde*M
     rm("M")
     # Garbage collection
     gc(verbose = FALSE)
   }

   if(Vs=="LoNL"){ 

     if(verbose){
       message("Computing N...")
     }

     # Calculate the summation in two steps to save memory: (1) Write down the first term, (2) add the remaining terms only if ngen2 - ngen1 >1

     if(ngen1!=0){  
       N = {((1-1/(2*nind*Ne_factor))^(ngen1))*nR^(ngen1)}
     }else{
       if(ngen2>1){
         N = matrix(0, ncol(nR), ncol(nR))
       }
     }

     # Only perform further summations if (ngen2-ngen1 > 1)
     
     if((ngen2 - ngen1) > 1){
     
       for (x in (ngen1 + 2):(ngen2)){
         N = N + {((1-1/(2*nind*Ne_factor))^(x-1))*nR^(x-1)}
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

    sdLoM<-RSpectra::eigs(t(UL)%*%Drift%*%UL, retain)
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
  
  

  if(is.na(palpha)){

    if(verbose){
      message("Estimating palpha...")
    }

    palpha<-optim(0, fit.model, balpha=balpha, LDalpha = LDalpha, nsnps=nsnps, UL=UL, DL=DL, L=L, ngen2=ngen2, ngen1=ngen1, nind=nind, tprojp=tprojp, pbar0=pbar0, pbar1=pbar1, pbar2=pbar2, nrep=nrep, LLonly=TRUE, Selec=Selec, Ne_factor=Ne_factor, verbose=verbose, method = "L-BFGS-B", lower = -2, upper =2, control = list(fnscale=-1, factr = 1e+11), hessian=TRUE)

    palpha_var<--1/palpha$hessian
    palpha<-palpha$par

  }else{
    palpha_var<-0
  }

  if(verbose){
    message("Fitting the final model...")
  }

  output<-fit.model(palpha=palpha, balpha=balpha, LDalpha = LDalpha, nsnps=nsnps, UL=UL, DL=DL, L=L, ngen2=ngen2, ngen1=ngen1, nind=nind, tprojp=tprojp, pbar0=pbar0, pbar1=pbar1, pbar2=pbar2, nrep=nrep, LLonly=FALSE, Selec=Selec, Ne_factor=Ne_factor, verbose=verbose)

  if(verbose){
    message("Calculating the estimate of Vw...")
  }

  if(method=="REML"){
    
    sigma2alpha<-summary(output$model)$varcomp[1,1]
    S<-matrix(0,2,2)

    if(is.na(balpha[1]) & is.na(balpha[2])){
        balpha<-summary(output$model, coef=TRUE)$coef.fixed[c("int", "pmq"),1]
        S<-output$model$Cfixed[c("int", "pmq"),c("int", "pmq")]
    }
    if(is.na(balpha[1]) & !is.na(balpha[2])){
       balpha[1]<-summary(output$model, coef=TRUE)$coef.fixed["int",1]
       S["int","int"]<-output$model$Cfixed
    }
    if(!is.na(balpha[1]) & is.na(balpha[2])){
        balpha[2]<-summary(output$model, coef=TRUE)$coef.fixed["pmq",1]
        S["pmq","pmq"]<-output$model$Cfixed
    }
  }
  if(method=="MCMC"){
    sigma2alpha<-posterior.mode(output$model$VCV[,1])
    balpha<-colMeans(output$model$Sol[,c("int", "pmq")])
  }

  X<-cbind(rep(1, length(pbar0)), 2*pbar0-1)

  if(LDalpha){
    TrV<-sum(DL^(2*(output$palpha+1)))*sigma2alpha
    aLa<-t(X%*%balpha)%*%L%*%X%*%balpha-sum(diag(t(X)%*%L%*%X%*%S))
  }else{
    TrV<-sum(diag(L)^(output$palpha+1))*sigma2alpha
    aLa<-t(X%*%balpha)%*%L%*%X%*%balpha-sum(diag(t(X)%*%L%*%X%*%S))
  }

  Vw_est<-TrV+aLa

  return(list(Vw_est=Vw_est, data=output$data, model=output$model, SC=output$SC, palpha=output$palpha, balpha=balpha, palpha_var=palpha_var, balpha_var=S, tprojp=if(save_tprojp){tprojp}else{NULL}, S=S, X=X, DL=ifelse(exists("DL"), DL, NA)))

}

logit<-function(x){
  log(x)-log1p(-x)
}

pmq_trans<-function(pmq, balpha_0){
  x<-sign(pmq)*abs(pmq)^balpha_0
  logit((x+1)/2)
}


alpha_distribution<-function(alpha, p, tprojp=NULL, logit=FALSE, save_model=FALSE){

  pmq<-2*p-1
  pq<-p*(1-p)/2

  palpha_fit<-function(palpha, alpha, pmq, pq, tprojp){

       if(!is.null(tprojp)){
          alpha<-t(alpha%*%tprojp)
          pmq<-t(pmq%*%tprojp)
          pq<-diag(t(tprojp)%*%(tprojp*pq^palpha))
       }else{
          pq<-pq^palpha
       }

       -logLik(lm(alpha~pmq, weights=1/pq))
  }

  logit_fit<-function(par, alpha, pmq, pq,  tprojp){

        balpha_0<-par[1]
        palpha<-par[2]

        pmq<-pmq_trans(pmq, balpha_0)

        if(!is.null(tprojp)){
          alpha<-t(alpha%*%tprojp)
          pmq<-t(pmq%*%tprojp)
          pq<-diag(t(tprojp)%*%(tprojp*pq^palpha))
        }else{
          pq<-pq^palpha
        }
        -logLik(lm(alpha~pmq-1, weights=1/pq))
  }

  if(logit){
    palpha<-optim(c(1, 1), logit_fit, alpha=alpha, pmq=pmq, pq=pq, tprojp=tprojp)$par

    balpha_0<-palpha[1]
    palpha<-palpha[2]
  
    pmq<-pmq_trans(pmq, balpha_0)

    if(!is.null(tprojp)){
      alpha<-t(alpha%*%tprojp)
      pmq<-t(pmq%*%tprojp)
      pq<-diag(t(tprojp)%*%(tprojp*pq^palpha))
    }else{
      pq<-pq^palpha
    }

    m1<-lm(alpha~pmq-1, weights=1/pq)
  }else{
    palpha<-optim(0, palpha_fit, alpha=alpha, pmq=pmq, pq=pq, tprojp=tprojp, method="Brent", lower=-5, upper=5)$par

    if(!is.null(tprojp)){
      alpha<-t(alpha%*%tprojp)
      pmq<-t(pmq%*%tprojp)
      pq<-diag(t(tprojp)%*%(tprojp*pq^palpha))
    }else{
      pq<-pq^palpha
    }

    m1<-lm(alpha~pmq, weights=1/pq)
  }  

  
  result<-list(palpha=palpha, 
               balpha_0=if(logit){balpha_0}else{coef(m1)["(Intercept)"]}, 
               balpha_1=coef(m1)["pmq"],
               sigma2alpha=summary(m1)$sigma^2,
               logit=logit,
               model=if(save_model){m1}else{NULL},
               X=model.matrix(m1),
               S=summary(m1)$cov.unscaled*summary(m1)$sigma^2)

  return(result)
}

plot_alpha_distribution<-function(alpha=NULL, p=NULL, tprojp=NULL, parameters, ...){

  pmq<-2*p-1
  pq<-p*(1-p)/2

  if(parameters$logit){
    pred<-pmq_trans(pmq, parameters$balpha_0)
    intercept<-0
  }else{
    pred<-pmq
    intercept<-parameters$balpha_0
  }
  slope<-parameters$balpha_1

  palpha<-parameters$palpha
  var<-pq^palpha

  if(is.null(tprojp)){

    pmq.x<-seq(min(pmq),max(pmq), length=1000)
    pq.x = (pmq.x+1)*(1-pmq.x)/8
    sd.x<-parameters$sigma2alpha*sqrt(pq.x^palpha)

    plot(alpha~pmq, ylab="alpha", xlab="p-q", ...)

    if(parameters$logit){
      pred.x<-pmq_trans(pmq.x, parameters$balpha_0)
    }else{
      pred.x<-pmq.x
    } 

    m.x<-intercept+slope*pred.x
    u.x<-qnorm(0.975, intercept+slope*pred.x, sd.x)
    l.x<-qnorm(0.025, intercept+slope*pred.x, sd.x)

    lines(m.x~pmq.x, col="red")
    lines(u.x~pmq.x, col="red", lty=2)
    lines(l.x~pmq.x, col="red", lty=2)
  
  }else{  

    alpha<-t(alpha%*%tprojp)
    pred<-t(pred%*%tprojp)
    var<-diag(t(tprojp)%*%(tprojp*pq^palpha))

    plot(alpha~pred, ylab="projected alpha", xlab="projected predictor", ...)

    pred.x<-seq(min(pred),max(pred), length=1000)

    m.x<-intercept+slope*pred.x
    lines(m.x~pred.x, col="red")

  }
}  


summaryG<-function(x){

  if(ncol(x)!=2 | any(!x%in%c(1:4))){
    stop("x should be an nx2 matrix with integers 1:4 corresponding to haplotypes AB (1), Ab (2), aB (3) and ab (4)")
  }

  n_hap<-length(x)

  pA<-sum(x==1 | x==2)/n_hap
  pB<-sum(x==1 | x==3)/n_hap
  # frequency of A and B alleles

  Lgp<-sum(x==1)/n_hap
  # frequency of AB haplotype

  Lgp<-Lgp-pA*pB
  # gametic phase LD

  Lngp<-sum((x[,1]==1 | x[,1]==2) & (x[,2]==1 | x[,2]==3))
  # number of individuals whose first haplotype contains an A and the second a B
  Lngp<-Lngp+sum((x[,2]==1 | x[,2]==2) & (x[,1]==1 | x[,1]==3))
  # add number of individuals whose first haplotype contains a B and the second an A
  Lngp<-Lngp/n_hap
  # frequency of individuals with A and B alleles on different haplotypes

  Lngp<-Lngp-pA*pB
  # non-gametic phase LD
 
  return(list(pA=pA, pB=pB, Lgp=Lgp, Lngp=Lngp))
}

sample_gamete<-function(xi, r){
   
   r_event<-rbinom(1, size=1, prob=r)

   if(r_event & all(c(1,4)%in%xi)){xi<-c(2,3)}
   if(r_event & all(c(2,3)%in%xi)){xi<-c(1,4)}

   h_event<-rbinom(1, size=1, prob=1/2)+1

   return(xi[h_event])
}

prob_novelA<-function(pA, pB, rD, nH){ 
  # pA and pB allele frequencies of reference alleles at each locus
  # rD measured as a correlation. 
  # nH the total number of parental haplotypes

  pa<-1-pA
  pb<-1-pB
  D<-rD*sqrt(pA*pa*pB*pb)


  P<-pA*pB+D
  Q<-pA*pb-D
  R<-pa*pB-D
  S<-pa*pb+D
  (P*S*(1-Q)^nH+P*S*(1-R)^nH+Q*R*(1-P)^nH+Q*R*(1-S)^nH)
}

prob_novelB<-function(P,Q,R,S, nH){ 
  # P Q R & S are the four haplotype frequencies
  # nH the total number of parental haplotypes


  P*S*(1-Q)^nH+P*S*(1-R)^nH+Q*R*(1-P)^nH+Q*R*(1-S)^nH
}

Q<-function(pbar1, pbar2){

   if(any(pbar2==0) | any(pbar2==1)){
     pbar2[which(pbar2==0 | pbar2==1)]<-NA
   }
   if(any(pbar1==0) | any(pbar1==1)){
     pbar1[which(pbar1==0 | pbar1==1)]<-NA
   }
   Q<-cov(t(pbar2-pbar1), use="pairwise.complete.obs")
   d<-rowMeans(pbar1*(1-pbar1), na.rm=TRUE)
   Q<-Q/outer(d,d)
   return(Q)
}

Sigma<-function(L,nR, N=Inf, nrep){

  n<-nrow(L)

  s<-(sum(abs(cov2cor(L))*nR)-n)/(n*(n-1))

  Sigma<-matrix(s, nrep, nrep)
  diag(Sigma)<-1+1/N

  return(Sigma)

}

est_Va_bc<-function(pbar1, pbar2, L, nR, nrep){

  Q<-Q(pbar1, pbar2)
  Sigma<-Sigma(L, nR, nrep=nrep)

  q<-Q[upper.tri(Q, diag=TRUE)]
  a<-Sigma[upper.tri(Sigma, diag=TRUE)]
  b<-diag(nrow(Q))[upper.tri(Q, diag=TRUE)]

  coef(lm(q~a+b-1))["a"]

}

extract_slim_data = function(Set_ID,                # The unique ID of the set of simulations that are controlled by a single R script
                             sim = 1,               # Each set can have multiple sims, but - on the cluster sim must always 1
                             unzip = FALSE,          # Should the SLiM output file be unzipped, read, and then zipped back?
                             slim_output_path,      # The directory where the SLiM outputs (for parents and experimental replicates) are stored (as .txt files)
                             sim_param_path,        # The path to the directory where the .csv file containing simulation parameters is stored
                             extract_genomes_path,  # The path to the python script that extracts genomes and mutations from SLim outputs (3_Extract_genomes.py)
                             extract_mut_path,      # The path to the python script that extracts mutations from SLim outputs (2_Extract_mulations.py)
                             mutations_path,        # The directory where extracted mutations are to be stored (temp files)
                             c_matrix_path,         # The directory where extracted genomes are to be stored (temp files)
                             n_sample=NULL,         # Number of individuals sampled from the parents' generation (useful if n_ind_exp is large)
                             randomise = TRUE,      # Optionally the reference allele can be randomised
                             delete_temp_files = TRUE, 
                             verbose=TRUE){     
  
  ################################
  ## Read simulation parameters ##
  ################################
      
  sim_params = read.csv(paste(sim_param_path, "/", Set_ID, "_Data.csv", sep = ""), header = T)
  # Restrict the data to current sim if the data frame has data from multiple sims
  if("sim"%in%colnames(sim_params)){sim_params = sim_params[sim_params$sim==sim,]}
  
  if(nrow(sim_params)==0){stop("No data found!")}
  
  # Load sim parameters
  
  n_ind_exp = sim_params$n_ind_exp     # Number of individuals in the experiment
  n_cages = sim_params$n_cages         # Number of replicates in the experiment
  if("ngen1"%in%colnames(sim_params)){ngen1 = sim_params$ngen1}else{ngen1 = 1}
  if("ngen2"%in%colnames(sim_params)){ngen2 = sim_params$ngen2}else{ngen2 = sim_params$ngen_expt + 1}
  ngen_expt = ngen2 - ngen1     # Number of generations over which allele frequency changes are calculated
  end_gen = sim_params$end_gen         # The generation number of the parents' generation

  
  # If n_sample is not provided extract the genomes of all individuals in the parents' generation
  
  if(is.null(n_sample)){n_sample = n_ind_exp}
  if(n_sample>n_ind_exp){stop("n_sample cannot be greater than n_ind_exp")}
  
  if(verbose){
    message("Reading the state of the population in the parent's generation...")
    message("Extracting mutations and genomes from the parents' generation...")
  }
  
  # Unzip the SLiM output file
  
  if(unzip){
    if(verbose)message("Unzipping the SLiM output file...")
    system(paste("gunzip", paste(slim_output_path, "/", Set_ID, "_sim", sim, "_output_parents.txt.gz", sep = "")))
  }

  system(paste("python", 
               extract_genomes_path,                        # Path of the python script (3_Extract_genomes.py)
               paste(slim_output_path, "/", Set_ID, "_sim", sim, "_output_parents.txt", sep = ""),  # Path of the .txt file containing the SLiM output for the parent's generation 
               paste(mutations_path, "/", Set_ID, "_sim", sim, "_mutations_parents.txt", sep = ""), # Path of the .txt output file containing the mutations in the parents' generation 
               paste(c_matrix_path,"/", Set_ID, "_sim", sim, "_c_matrix_parents.csv", sep = ""),    # Path of the .csv output file containing the c matrix for genomes in the parents' generation
               n_sample))                                                                 # Number of individuals to be sampled randomly to construct the c matrix (just for space issues in case n_ind_exp is very large). Typically should be set to same as n_ind_exp
  
  # Rezip the SLiM output file
  
  if(unzip){
    if(verbose)message("Re-zipping the unzipped SLiM output file...")
    system(paste("gzip", paste(slim_output_path, "/", Set_ID, "_sim", sim, "_output_parents.txt", sep = "")))
  }

  
  # Read genomes
  c_genome = read.csv(paste(c_matrix_path,  "/", Set_ID, "_sim", sim, "_c_matrix_parents.csv", sep =""), header=F) # as.integer done to avoid scientific notation
  
  # Delete the .csv file
  if(delete_temp_files){system(paste("rm", paste(c_matrix_path,  "/", Set_ID, "_sim", sim, "_c_matrix_parents.csv", sep =""), sep = " "))}
  
  c_genome = as.matrix(c_genome)
  
  # Convert genome data into individual data (rows are individuals and columns are allele counts {0,1 or 2} at various sites). 
  # Note that genome 1 and genome 2 are from individual 1; 3 and 4 are from individual 2, and so on
  
  n0_individuals = nrow(c_genome)/2

  # If one samples individuals from the parents' generation while building the c matrix (i.e. when sample_size is less than n_ind_exp), the sample may not contain some low frequency mutations, i.e. some loci are not segregating in the sample, but are in the parents' population
  
  # identify the loci that are missing in the sample
  
  retained_loci = which(colSums(c_genome)!=0)
  
  # Trim the c matrix to contain only the retained loci, i.e. the loci that are segregating in the sample
  
  c_genome = c_genome[,retained_loci]
  n_sites = ncol(c_genome)
    
  ## Read the list of mutations output generated by SLiM and subsequently cropped out as a separate file by Python (for the parents' generation)

  pbar0 = colMeans(c_genome)
  # vector of starting allele frequencies

  mutations_0 = read.table(paste(mutations_path, "/", Set_ID, "_sim", sim, "_mutations_parents.txt", sep = ""), header=T, sep = " ")
  
  # Delete the .txt file
  if(delete_temp_files){system(paste("rm ", mutations_path, "/", Set_ID, "_sim", sim, "_mutations_parents.txt", sep = ""))}
  
  # Sort mutations based on the Temp_ID, so that the order of loci matches the order in c_ind
  
  mutations_0 = mutations_0[order(mutations_0$Temp_ID),]
  
  # Trim mutations to contain only the retained loci, i.e. the loci that are segregating in the sample
  
  mutations_0 = mutations_0[retained_loci,]
  
  list_alpha = 2*(mutations_0$s)   # Vector of alphas
  SNPs = mutations_0$Position      # Vector of positions of mutations in the parents
  
  ############################################################################
  ############################################################################
  ####### Calculate allele frequencies in the experiment for each cage #######
  ############################################################################
  ############################################################################  
  
  # Create empty vector to create the data frame containing the following variables as columns:
  # 1. Raw delta P
  # 2. Projected delta P
  # 3. Cage ID (replicate)
  # 4. Locus ID
  
  d_proj = c()
  d_raw = c()
  P_matrix = c()  
  
  for (cage in (1:n_cages)){
    
    ###############################
    ###### extract mutations ######
    ###############################
    
    
    # There are two command line arguments (1. Path of the SLiM output file, 2. Path where mutations are to be written) 
    
    if(verbose){
      message(paste("Extracting mutations and storing allele frequencies in cage ", cage, " of simulation ", sim, "...", sep = ""))
    }

    for (gen in (end_gen + 1):(end_gen + ngen2)){
      
      # Unzip SLiM output files
      
      if(unzip){
        if(verbose){message("Unzipping the SLiM output file...")}
        system(paste("gunzip", paste(slim_output_path, "/", Set_ID, "_sim", sim, "_cage", cage, "_output_experiment_", as.integer(gen), ".txt", sep = "")))
      }
      
      system(paste("python", 
                   extract_mut_path, 
                   paste(slim_output_path, "/", Set_ID, "_sim", sim, "_cage", cage, "_output_experiment_", as.integer(gen), ".txt", sep = ""), 
                   paste(mutations_path, "/", Set_ID, "_sim", sim, "_cage", cage, "_mutations_", as.integer(gen), ".txt", sep = ""))) # as.integer done to avoid scientific notation
      
      # Re-zip unzipped SLiM files
      if(unzip){
        if(verbose){message("Re-zipping the unzipped SLiM output file...")}
        system(paste("gzip", paste(slim_output_path, "/", Set_ID, "_sim", sim, "_cage", cage, "_output_experiment_", as.integer(gen), ".txt", sep = "")))
      }
    }
    
    
    ### Create an empty matrix to store allelic frequencies in each generation
    
    # If some of the loci get fixed/lost in subsequent generations, NAs should be inserted
    
    # Create an empty vector to store allele frequencies
    
    P = c()
    
    # store the frequencies in the parent's generation in P
    # Frequency = (Number of genomes)/(2*popsize)
    
    P = cbind(P, mutations_0$Number/(2*n_ind_exp))
    
    ### Loop through the generations of the experiment identifying mutations that were present in the parents' generation (using permanent IDs) and recording their frequencies
    
    for (gen in (end_gen+1):(end_gen+ngen2)){
      
      # Read the file storing mutation information
      mut = read.csv(paste(mutations_path, "/", Set_ID, "_sim", sim, "_cage", cage, "_mutations_", as.integer(gen), ".txt", sep = ""), sep = " ") # as.integer gets rid of scientific notation
      
      # Delete the .txt file
      if(delete_temp_files){system(paste("rm ", mutations_path, "/", Set_ID, "_sim", sim, "_cage", cage, "_mutations_", as.integer(gen), ".txt", sep = ""))}
      
      # Sort mutations based on the Temp_ID, so that the order of loci matches the order in the C and L matrices
      
      mut = mut[order(mut$Temp_ID),]
      
      # Create an empty vector to store frequencies of mutations in the current generation
      freq = c()
      
      # Loop through the permanent IDs of  mutations segregating in end_gen (parents' generation)
      # i.e. Loop through Permanent IDs in mutations_0
      # Check if each mutation is present in the current generation
      # If present, record the frequency in freq, otherwise add either 0 or 1 to freq using the round() function
      
      for(mutation in mutations_0$Permanent_ID){
        if(mutation %in% mut$Permanent_ID){freq = c(freq, (mut[which(mut$Permanent_ID==mutation),]$Number)/(2*n_ind_exp))
        }else{
          freq = c(freq, round(P[which(mutations_0$Permanent_ID==mutation), gen - end_gen]))
          
        }
      }
      
      # Add the vector freq to P
      
      P = cbind(P, freq)
      
    }
   
    P_matrix = rbind(P_matrix, P)
    
  }
  
  
  ##### Calculate matrices of allele frequencies with rows as replicates ######
  
  pbar1 = matrix(NA, nrow = n_cages, ncol = n_sites)
  pbar2 = matrix(NA, nrow = n_cages, ncol = n_sites)
  for (i in 1:n_cages){
    pbar1[i,] = P_matrix[((i-1)*n_sites + 1):((i-1)*n_sites + n_sites), (ngen1 + 1)] ## Matrix of frequencies at the start of the experiment (ie F1 generation)
  }
  
  for (i in 1:n_cages){
    pbar2[i,] = P_matrix[((i-1)*n_sites + 1):((i-1)*n_sites + n_sites), (ngen2 + 1)] ## Matrix of frequencies at the end of the experiment
    
    
  }
    
  if(randomise){
    
      #########################################################
      ######## Randomise the reference allele in c_ind ########
      #########################################################
      
      # Randomly change the reference allele
      # This can be done as follows:
      # Add 0s to the allele counts of those alleles that stay the same, and -2 to those alleles that are to be switched
      # Then take a mod
      # Remember c_ind is a matrix of 0s 1s and 2s
      
      if(verbose){
        message("Randomising reference alleles in c_ind...")
      }
        
      # Generate a random vector of 0s (for no change) and -1s (for loci where the reference allele is to be switched)

      ran_vect = sample(c(0, -1), ncol(c_genome),  replace = T) 
      
      # Create a matrix with with as many rows as c_ind. 
      # Each row of this matrix should be made up of two times ran_vect (since we are working with allele counts, not frequencies). 
      # Because the same changes need to be applied to each individual
 
      ran_matrix = t(matrix(ran_vect, nrow = ncol(c_genome), ncol = nrow(c_genome)))
      
      # Calculate the allele counts of the new (randomised) reference alleles
      
      c_genome = abs(c_genome + ran_matrix)

      ##############################################################
      ######## Randomise the reference allele in list_alpha ########
      ##############################################################
      
      if(verbose){
        message("Randomising reference alleles in list_alpha...")
      }

      # Create a new vector of -1s (wherever ran_vect has -1) and 1s (wherever ran_vect has 0)
      ran_vect_alpha = ifelse(ran_vect==-1, -1, 1)
      
      list_alpha = list_alpha*ran_vect_alpha
      
      
      ################################################################
      ######## Randomise the reference allele in pbar1 and pbar2 #####
      ################################################################
      
      # Randomly change the reference allele
      # This can be done as follows:
      # Add 0s to the frequencies of those alleles that stay the same, and -1 to those alleles that are to be switched
      # Then take a mod
      
      if(verbose){
        message("Randomising reference alleles in pbar1 and pbar2...")
      }

      # Create a matrix with with as many rows as pbar1. Each row of this matrix should be made up of ran_vect (computed while randomising c_ind). 
      # Because the same changes need to be applied to each cage
      
      ran_matrix_pbar = t(matrix(ran_vect, nrow = ncol(pbar1), ncol = nrow(pbar1)))
      
      # Calculate the frequencies of the new (randomised) reference alleles

      pbar0 = abs(pbar0 + ran_vect)     
      pbar1 = abs(pbar1 + ran_matrix_pbar)
      pbar2 = abs(pbar2 + ran_matrix_pbar)      
    
  }

  return(list(c_genome=c_genome, list_alpha=list_alpha, SNPs=SNPs, ngen1=ngen1, ngen2=ngen2,  pbar0=pbar0, pbar1=pbar1, pbar2=pbar2, sim_params=sim_params[,1:20]))
  
}

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
  
  return(list(L=L, Ltilde=Ltilde, nR=nR, svdL=if(compute_svdL | LDalpha){list(UL=UL, DL=DL)}else{NULL}, seg_sites=seg_sites, seg_sites_neu=seg_sites_neu, seg_sites_ben=seg_sites_ben, seg_sites_del=seg_sites_del, mean_diversity=mean_diversity, va_true=va_true, vA_true=vA_true, vA_alpha_emp=vA_alpha_emp, parameters=alpha_properties[1:5]))

}

analyse_sim = function(Set_ID,                # The unique ID of the set of simulations that are controlled by a single R script
                       sim = 1,               # Each set can have multiple sims, but - on the cluster sim must always 1
                       unzip = FALSE,          # Should the SLiM output file be unzipped, read, and then zipped back?
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
                       proj = "BLoM",           # projection type for allele frequencies: "LoM", "BLoM", "L" or "N"
                       LDalpha = FALSE,       # Should L or diag(L) be considered while modelling distribution of alphas
                       pa = 1,
                       Vs = "LoNL",           # "L" or "LoNL"
                       method="REML",         # Can be "REML" or "MCMC"
                       palpha = NA,           # If NA palpha is estimated using optim()
                       balpha = c(NA, NA),    # If c(NA,NA) both bedelta intercept and slope are estimated
                       AtleastOneRecomb=FALSE, 
                       all.gp = FALSE,        # Ltilde = L'+L''(r/(1-r)) if all.gp=T L'' is assumed 0 and L'=L. 
                       verbose = TRUE
                       ){
  
  ### Extract sim data ###
  
  if(verbose){message("Extracting simulation data...")}

  sim_data = extract_slim_data(Set_ID = Set_ID,
                               sim = sim,
                               unzip = unzip,
                               slim_output_path = slim_output_path, 
                               sim_param_path = sim_param_path,
                               extract_genomes_path = extract_genomes_path, 
                               extract_mut_path = extract_mut_path,
                               mutations_path = mutations_path, 
                               c_matrix_path = c_matrix_path, 
                               randomise = randomise)
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
               tol = sqrt(.Machine$double.eps))


  vA_est = m1$Vw_est 
  palpha_est = m1$palpha 
  palpha_var_est = m1$palpha_var
  balpha_intercept_est = m1$balpha[1]
  balpha_slope_est = m1$balpha[2]
  balpha_var_est = paste(m1$balpha_var[1,1], m1$balpha_var[2,2], m1$balpha_var[1,2], sep = "_")
  sigma2alpha_est = summary(m1$model)$varcomp[1,1]
    
  ### Save file ###
  
  if(verbose){message("Saving data...")}
 
  sim_params = sim_data$sim_params

  analysis_data = data.frame("proj"=proj, "LDalpha"=LDalpha, "pa"=pa, "Vs"=Vs, "randomise"=randomise, "palpha_method"=palpha, "balpha_method"=paste(balpha[1], balpha[2], sep="_"), "va_true"=parents_info$va_true, "vA_true"=parents_info$vA_true, "vA_est"=vA_est, "vA_alpha_emp"=parents_info$vA_alpha_emp, "palpha_emp"=parents_info$parameters$palpha, "balpha_intercept_emp"=parents_info$parameters$balpha_0, "balpha_slope_emp"=parents_info$parameters$balpha_1, "sigma2alpha_emp"=parents_info$parameters$sigma2alpha, "palpha_est"=palpha_est, "palpha_var_est"=palpha_var_est, "balpha_intercept_est"=balpha_intercept_est, "balpha_slope_est"=balpha_slope_est, "balpha_var_est"=balpha_var_est, "sigma2alpha_est"=sigma2alpha_est, "seg_sites"=parents_info$seg_sites, "seg_sites_neu"=parents_info$seg_sites_neu, "seg_sites_ben"=parents_info$seg_sites_ben, "seg_sites_del"=parents_info$seg_sites_del, "mean_diversity"=parents_info$mean_diversity, "all.gp" = all.gp)

  analysis_data = cbind(sim_params, analysis_data)
  write.table(rbind(names(analysis_data), analysis_data), file = paste(output_path, "/", Set_ID, "_sim_", sim, "_Data_analysis.csv", sep = ""),col.names = FALSE, row.names = FALSE, sep = ",")
  
  return(analysis_data)
  
}
