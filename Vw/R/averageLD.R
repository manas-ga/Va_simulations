
averageLD<-function(L, nR, nrep){
  
  # Equation 55 in Buffalo & Coop
  n<-nrow(L)
  
  s<-(sum((cov2cor(L)^2)*nR)-n)/(n*(n-1))
  # Equation 55 Buffalo & Coop
  
  A<-matrix(s, nrep, nrep)
  
  return(A)
  
}
