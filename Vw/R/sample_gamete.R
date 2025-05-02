sample_gamete<-function(xi, r, n=1){
  
  r_event<-rbinom(n, size=1, prob=r)
  h_event<-rbinom(n, size=1, prob=1/2)+1
  
  gametes<-1:n

  gametes<-xi[h_event]

  gametes[which(r_event & all(c(1,4)%in%xi) & h_event==1)]<-2
  gametes[which(r_event & all(c(1,4)%in%xi) & h_event==2)]<-3
  gametes[which(r_event & all(c(2,3)%in%xi) & h_event==1)]<-1
  gametes[which(r_event & all(c(2,3)%in%xi) & h_event==2)]<-4
  
  return(gametes)

}

