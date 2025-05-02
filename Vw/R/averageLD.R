
averageLD<-function(L, nR, nrep, selected=NULL, ngen=1){
  
  # Equation 55 in Buffalo & Coop

  if(any(dim(L)!=dim(nR))){
    stop("L and nR should have the same dimensions")
  }

  n<-nrow(L)
  
  LD<-(cov2cor(L)^2)*nR^(2*ngen)

  if(is.null(selected)){
    s<-(sum(LD)-n)/(n*(n-1))
  }else{
    if(min(selected)<1 | max(selected)>nrow(L)){
      stop("selected should be loci between 1 and the dimension of L")
    }
    s<-mean(LD[selected, -selected])
  }  

  A<-matrix(s, nrep, nrep)
  
  return(A)
  
}
