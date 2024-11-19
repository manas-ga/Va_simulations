
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
