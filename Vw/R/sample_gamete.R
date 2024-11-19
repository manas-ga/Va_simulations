sample_gamete<-function(xi, r){
  
  r_event<-rbinom(1, size=1, prob=r)
  
  if(r_event & all(c(1,4)%in%xi)){xi<-c(2,3)}
  if(r_event & all(c(2,3)%in%xi)){xi<-c(1,4)}
  
  h_event<-rbinom(1, size=1, prob=1/2)+1
  
  return(xi[h_event])
}