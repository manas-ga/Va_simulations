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
