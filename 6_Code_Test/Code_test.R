rm(list = ls())



### Load my simulation data
load("C:/Users/msamant/Downloads/test3.RData")

##### Calculate matrices of allele frequencies with rows as replicates ######

pbar1 = matrix(NA, nrow = n_cages, ncol = n_sites)
pbar2 = matrix(NA, nrow = n_cages, ncol = n_sites)
for (i in 1:n_cages){
  pbar1[i,] = P_matrix[((i-1)*n_sites + 1):((i-1)*n_sites + n_sites),2]
}

for (i in 1:n_cages){
  pbar2[i,] = P_matrix[((i-1)*n_sites + 1):((i-1)*n_sites + n_sites),5]
  
  
}


#### Define the function




Vw_model<-function(C0,          # parental genotypes (rows individuals, columns loci, coded as 0, 1/2 or 1) 
                   nR,          # matrix of non-recombinant probabilities between loci
                   pbar1,       # vector of allele frequencies at time-point 1
                   ngen1=1,     # number of generations between parents and time-point 1
                   pbar2,       # vector of allele frequencies at time-point 2
                   ngen2,       # number of generations between parents and time-point 2
                   nind,        # population size in each replicate
                   proj="BLoM", # projection type for allele frequencies: "LoM", "BLoM", "L" or "N"
                   LDdelta,
                   pa,
                   pdelta=0,
                   Vs,
                   method="REML",
                   L=NULL,    # list with elements UL and DL
                   svdL=NULL,    # list with elements UL and DL
                   tol=sqrt(.Machine$double.eps),
                   nrep=10
){
  
            if(is.null(L)){
              L<-cov(C0)*(nind-1)/(nind)
            }
            
            if(!proj%in%c("LoM", "BLoM", "L", "N")){stop("proj must be one of 'LoM', 'L', 'N'")}
            if(!Vs%in%c("LoNL", "L")){stop("Vs must be either 'LoNL' or 'L'")}
            if(!method%in%c("REML", "MCMC")){stop("method must be either 'REML' or 'MCMC'")}
            
            #################################
            # calculate projection matrices #
            #################################
            nsnps<-ncol(C0)
            
            if(proj=="L" | proj=="BLoM" | LDdelta){
              
              if(is.null(svdL)){
                svdC<-svd(t(scale(sqrt(1/nind)*C0, scale=FALSE)))
                retain<-sum(svdC$d>tol)
                U<-svdC$u[,1:retain]
                D<-diag(svdC$d[1:retain])
              }else{
                U<-svdL$UL
                D<-diag(svdL$DL)
                retain<-ncol(U)
              }
              if(LDdelta){
                DL<-D
              }
            }
            
            if(proj=="LoM" | proj=="BLoM"){ 
              
              M<-Reduce('+', sapply(1+ngen1:(ngen2-ngen1), function(x){((1-1/(2*nind))^(x-1))*(1/nind)*nR^x}, simplify=FALSE))
              
              if(proj=="LoM"){
                sdLoM<-RSpectra::eigs(L*M, min(nind, nsnps))
                retain<-sum(sqrt(sdLoM$values)>tol)
                U<-sdLoM$vectors[,1:retain]
                D<-diag(sqrt(sdLoM$values[1:retain]))
              }
              
              if(proj=="BLoM"){ 
                sdLoM<-RSpectra::eigs(t(U)%*%(L*M)%*%U, ncol(U))
                retain<-sum(sqrt(sdLoM$values)>tol)
                U<-U%*%sdLoM$vectors[,1:retain]
                D<-diag(sqrt(sdLoM$values[1:retain]))
              }   
            } 
            
            
            if(proj=="N"){
              
              projp<-diag(nrow(pbar1)) #### changed from length(pbar1) to nrow(pbar1)
              
            }else{
              
              projp<-U%*%diag(diag(D)^(-pa))
              
            }
            
            if(LDdelta){
              if(pdelta==0){
                covp<-diag(nsnps)
              }else{
                covp<-UL%*%diag(diag(DL)^(2*pdelta))%*%t(UL)
              }  
            }else{
              covp<-diag(diag(L)^pdelta)
            }
            
            if(Vs=="LoNL"){
              N<-Reduce('+', sapply(1+ngen1:(ngen2-ngen1), function(x){((1-1/(2*nind))^(x-1))*nR^(x-1)}, simplify=FALSE))
              
              covp<-(L*N)%*%covp%*%(L*N)
            }
            if(Vs=="L"){
              covp<-(L*(ngen2 - ngen1))%*%covp%*%(L*(ngen2 - ngen1))
            }  
            
            
            SC<-t(projp)%*%covp%*%projp
            
            attr(SC, "INVERSE")<-FALSE
            dimnames(SC) <- list(1:nrow(SC),1:ncol(SC))  # used for full-form matrices, ## changing "retain" to "nrow(SC)" as retain does not exist when proj = "N" 
            
            pbar1_proj<-pbar1%*%projp 
            pbar2_proj<-pbar2%*%projp 
            
            dat.gaussian<-data.frame(delta=c(pbar2_proj-pbar1_proj), locus=gl(ncol(pbar1_proj),nrep,ncol(pbar1_proj)*nrep), rep=gl(nrep,1,ncol(pbar1_proj)*nrep))
            
            ##############
            # Fit models #
            ##############
            
            prior<-list(R=list(V=1, nu=0), G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)))
            
            if(method=="REML"){
              m1<-asreml(delta~1, random = ~vm(locus, SC), data=dat.gaussian)
            }
            if(method=="MCMC"){
              
              invSC<-solve(t(projp)%*%covp%*%projp)
              invSC <- as(invSC, "sparseMatrix") 
              attr(invSC, "rowNames") <- 1:retain
              attr(invSC, "colNames") <- 1:retain
              
              m1<-MCMCglmm(delta~1, random=~locus, data=dat.gaussian, ginverse=list(locus=invSC), family="gaussian", pr=TRUE, prior=prior)
            }
            
            if(LDdelta){
              TrV<-sum(diag(DL)^(2*(pdelta+1)))
            }else{
              TrV<-sum(diag(L%*%diag(diag(L)^pdelta)))
            }
            
            if(method=="REML"){
              Vw_est<-summary(m1)$varcomp[1,1]*TrV
            }
            if(method=="MCMC"){
              Vw_est<-posterior.mode(m1$VCV[,1])*TrV
            }
            
            return(list(Vw_est=Vw_est, data=dat.gaussian, model=m1, SC=SC))
            
}            





# Corrections in the function

# Line 72: changed from length(pbar1) to nrow(pbar1)
# Line 103: changed "retain" to "nrow(SC)" as retain does not exist when proj = "N"
# Line 105: changed pbar1 to t(pbar1) to make the multiplication conformable
# Line 106: changed pbar2 to t(pbar2) to make the multiplication conformable



Vw_model(C0 = c_ind_ret/2,          # parental genotypes (rows individuals, columns loci, coded as 0, 1/2 or 1) 
       nR = NRF,          # matrix of non-recombinant probabilities between loci
       pbar1 = pbar1,       # vector of allele frequencies at time-point 1
       ngen1=1,     # number of generations between parents and time-point 1
       pbar2 = pbar2,       # vector of allele frequencies at time-point 2
       ngen2 = 4,       # number of generations between parents and time-point 2
       nind = 1000,        # population size in each replicate
       proj="LoM", # projection type for allele frequencies: "LoM", "BLoM", "L" or "N"
       LDdelta = F,
       pa = 1,
       pdelta=0,
       Vs = "L",
       method="REML",
       L=NULL,    # list with elements UL and DL
       svdL=NULL,    # list with elements UL and DL
       tol=sqrt(.Machine$double.eps),
       nrep = 10)

###############################################################################################################
################################### Running things outside a function #########################################
###############################################################################################################


rm(list = ls())



### Load my simulation data
load("C:/Users/msamant/Downloads/test3.RData")

##### Calculate matrices of allele frequencies with rows as replicates ######

pbar1 = matrix(NA, nrow = n_cages, ncol = n_sites)
pbar2 = matrix(NA, nrow = n_cages, ncol = n_sites)
for (i in 1:n_cages){
  pbar1[i,] = P_matrix[((i-1)*n_sites + 1):((i-1)*n_sites + n_sites),2]
}

for (i in 1:n_cages){
  pbar2[i,] = P_matrix[((i-1)*n_sites + 1):((i-1)*n_sites + n_sites),5]
  
  
}


C0 = c_ind_ret/2          # parental genotypes (rows individuals, columns loci, coded as 0, 1/2 or 1) 
nR = NRF          # matrix of non-recombinant probabilities between loci
pbar1 = pbar1       # vector of allele frequencies at time-point 1
ngen1=1     # number of generations between parents and time-point 1
pbar2 = pbar2       # vector of allele frequencies at time-point 2
ngen2 = 4       # number of generations between parents and time-point 2
nind = 1000        # population size in each replicate
proj="LoM" # projection type for allele frequencies: "LoM", "BLoM", "L" or "N"
LDdelta = F
pa = 1
pdelta=0
Vs = "L"
method="REML"
L=NULL    # list with elements UL and DL
svdL=NULL    # list with elements UL and DL
tol=sqrt(.Machine$double.eps)
nrep = 10



if(is.null(L)){
  L<-cov(C0)*(nind-1)/(nind)
}

if(!proj%in%c("LoM", "BLoM", "L", "N")){stop("proj must be one of 'LoM', 'L', 'N'")}
if(!Vs%in%c("LoNL", "L")){stop("Vs must be either 'LoNL' or 'L'")}
if(!method%in%c("REML", "MCMC")){stop("method must be either 'REML' or 'MCMC'")}

#################################
# calculate projection matrices #
#################################
nsnps<-ncol(C0)

if(proj=="L" | proj=="BLoM" | LDdelta){
  
  if(is.null(svdL)){
    svdC<-svd(t(scale(sqrt(1/nind)*C0, scale=FALSE)))
    retain<-sum(svdC$d>tol)
    U<-svdC$u[,1:retain]
    D<-diag(svdC$d[1:retain])
  }else{
    U<-svdL$UL
    D<-diag(svdL$DL)
    retain<-ncol(U)
  }
  if(LDdelta){
    DL<-D
  }
}

if(proj=="LoM" | proj=="BLoM"){ 
  
  M<-Reduce('+', sapply(1+ngen1:(ngen2-ngen1), function(x){((1-1/(2*nind))^(x-1))*(1/nind)*nR^x}, simplify=FALSE))
  
  if(proj=="LoM"){
    sdLoM<-RSpectra::eigs(L*M, min(nind, nsnps))
    retain<-sum(sqrt(sdLoM$values)>tol)
    U<-sdLoM$vectors[,1:retain]
    D<-diag(sqrt(sdLoM$values[1:retain]))
  }
  
  if(proj=="BLoM"){ 
    sdLoM<-RSpectra::eigs(t(U)%*%(L*M)%*%U, ncol(U))
    retain<-sum(sqrt(sdLoM$values)>tol)
    U<-U%*%sdLoM$vectors[,1:retain]
    D<-diag(sqrt(sdLoM$values[1:retain]))
  }   
} 


if(proj=="N"){
  
  projp<-diag(nrow(pbar1)) #### changed from length(pbar1) to nrow(pbar1)
  
}else{
  
  projp<-U%*%diag(diag(D)^(-pa))
  
}

if(LDdelta){
  if(pdelta==0){
    covp<-diag(nsnps)
  }else{
    covp<-UL%*%diag(diag(DL)^(2*pdelta))%*%t(UL)
  }  
}else{
  covp<-diag(diag(L)^pdelta)
}

if(Vs=="LoNL"){
  N<-Reduce('+', sapply(1+ngen1:(ngen2-ngen1), function(x){((1-1/(2*nind))^(x-1))*nR^(x-1)}, simplify=FALSE))
  
  covp<-(L*N)%*%covp%*%(L*N)
}
if(Vs=="L"){
  covp<-(L*(ngen2 - ngen1))%*%covp%*%(L*(ngen2 - ngen1))
}  


SC<-t(projp)%*%covp%*%projp

attr(SC, "INVERSE")<-FALSE
dimnames(SC) <- list(1:nrow(SC),1:ncol(SC))  # used for full-form matrices, ## changing "retain" to "nrow(SC)" as retain does not exist when proj = "N" 

pbar1_proj<-pbar1%*%projp 
pbar2_proj<-pbar2%*%projp 

dat.gaussian<-data.frame(delta=c(pbar2_proj-pbar1_proj), locus=gl(ncol(pbar1_proj),nrep,ncol(pbar1_proj)*nrep), rep=gl(nrep,1,ncol(pbar1_proj)*nrep))

##############
# Fit models #
##############

prior<-list(R=list(V=1, nu=0), G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)))

if(method=="REML"){
  m1<-asreml(delta~1, random = ~vm(locus, SC), data=dat.gaussian)
}
if(method=="MCMC"){
  
  invSC<-solve(t(projp)%*%covp%*%projp)
  invSC <- as(invSC, "sparseMatrix") 
  attr(invSC, "rowNames") <- 1:retain
  attr(invSC, "colNames") <- 1:retain
  
  m1<-MCMCglmm(delta~1, random=~locus, data=dat.gaussian, ginverse=list(locus=invSC), family="gaussian", pr=TRUE, prior=prior)
}

if(LDdelta){
  TrV<-sum(diag(DL)^(2*(pdelta+1)))
}else{
  TrV<-sum(diag(L%*%diag(diag(L)^pdelta)))
}

if(method=="REML"){
  Vw_est<-summary(m1)$varcomp[1,1]*TrV
}
if(method=="MCMC"){
  Vw_est<-posterior.mode(m1$VCV[,1])*TrV
}

#return(list(Vw_est=Vw_est, data=dat.gaussian, model=m1, SC=SC))

Vw_est
