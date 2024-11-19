
prob_novelB<-function(P,Q,R,S, nH){ 
  # P Q R & S are the four haplotype frequencies
  # nH the total number of parental haplotypes
  
  
  P*S*(1-Q)^nH+P*S*(1-R)^nH+Q*R*(1-P)^nH+Q*R*(1-S)^nH
}