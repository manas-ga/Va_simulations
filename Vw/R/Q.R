Q<-function(pbar1, pbar2, exact=FALSE){
  
  # Equation 16 in Buffalo & Coop
  if(any(pbar2==0) | any(pbar2==1)){
    pbar2[which(pbar2==0 | pbar2==1)]<-NA
  }
  if(any(pbar1==0) | any(pbar1==1)){
    pbar1[which(pbar1==0 | pbar1==1)]<-NA
  }
  if(exact){
    p<-colMeans(sqrt(pbar1*(1-pbar1)))
    Q<-cov(t(pbar2)/p-t(pbar1)/p, use="pairwise.complete.obs")
  }else{
    d<-rowMeans(pbar1*(1-pbar1), na.rm=TRUE)  
    Q<-cov(t(pbar2-pbar1), use="pairwise.complete.obs")
    Q<-Q/(outer(d,d, "+")/2)
    # taken the average p*(1-p) across the pair
  }   
  return(Q)
}