
Sigma<-function(L, nR, nrep){
  
  # Equation 12 in Buffalo & Coop (multiplied by 2)
  n<-nrow(L)
  
  s<-(sum((cov2cor(L)^2)*nR)-n)/(n*(n-1))
  # Equation 55 Buffalo & Coop
  
  Sigma<-matrix(s, nrep, nrep)
  
  return(Sigma)
  
}
