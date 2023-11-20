rm(list = ls())

functions_only=TRUE
rmarkdown::render(file.path(if(Sys.info()["nodename"]=="sce-bio-c04553"){"~/Work/Va_simulations"}else{""} , "6_Code_Test/Vw.Rmd"))

### Load my simulation data
load("~/Downloads/test3.RData")

rm(list=ls()[-which(ls()%in%c("c_ind_ret", "NRF", "n_cages", "n_sites", "P_matrix", "run", "Vw_model"))])

##### Calculate matrices of allele frequencies with rows as replicates ######

pbar1 = matrix(NA, nrow = n_cages, ncol = n_sites)
pbar2 = matrix(NA, nrow = n_cages, ncol = n_sites)
for (i in 1:n_cages){
  pbar1[i,] = P_matrix[((i-1)*n_sites + 1):((i-1)*n_sites + n_sites),2]
}

for (i in 1:n_cages){
  pbar2[i,] = P_matrix[((i-1)*n_sites + 1):((i-1)*n_sites + n_sites),5]
  
  
}

nseq<-100

LL<-Vw_est<-1:nseq
pdelta<-seq(-2.5,-1.5, length=nseq)

for(i in 1:nseq){

  m1<-Vw_model(C0 = c_ind_ret/2,          # parental genotypes (rows individuals, columns loci, coded as 0, 1/2 or 1) 
         nR = NRF,          # matrix of non-recombinant probabilities between loci
         pbar1 = pbar1,       # vector of allele frequencies at time-point 1
         ngen1=1,     # number of generations between parents and time-point 1
         pbar2 = pbar2,       # vector of allele frequencies at time-point 2
         ngen2 = 4,       # number of generations between parents and time-point 2
         nind = 1000,        # population size in each replicate
         proj="BLoM", # projection type for allele frequencies: "LoM", "BLoM", "L" or "N"
         LDdelta = TRUE,
         pa = 1,
         pdelta=pdelta[i],
         Vs = "LoNL",
         method="REML",
         L=NULL,    # list with elements UL and DL
         svdL=NULL,    # list with elements UL and DL
         tol=sqrt(.Machine$double.eps))

  LL[i]<-m1$model$loglik
  Vw_est[i]<-m1$Vw_est
}

par(mfrow=c(2,1))
plot(LL~pdelta)
plot(Vw_est~pdelta)


