
Q<-function(pbar1, pbar2){
  
  # Equation 16 in Buffalo & Coop
  if(any(pbar2==0) | any(pbar2==1)){
    pbar2[which(pbar2==0 | pbar2==1)]<-NA
  }
  if(any(pbar1==0) | any(pbar1==1)){
    pbar1[which(pbar1==0 | pbar1==1)]<-NA
  }
  Q<-cov(t(pbar2-pbar1), use="pairwise.complete.obs")
  d<-rowMeans(pbar1*(1-pbar1), na.rm=TRUE)
  Q<-Q/(outer(d,d, "+")/2)
  # taken the average p*(1-p) across the pair 
  return(Q)
}