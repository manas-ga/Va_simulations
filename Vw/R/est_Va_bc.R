est_Va_bc<-function(pbar1, pbar2, L, nR, selected=NULL, exact=FALSE, ngen=1){
  
  # Equations 20 and 21 in Buffalo & Coop
  # should also scale a by the change in SSH
  
  nrep = nrow(pbar1)

  if(any(dim(pbar1)!=dim(pbar2))){
    stop("pbar1 and pbar2 should have the same dimension")
  }
  
  Q<-Q(pbar1, pbar2, exact=exact)

  if(any(dim(L)!=dim(nR))){
    stop("L and nR should have the same dimensions")
  }

  A<-averageLD(L, nR, nrep=nrep, selected=selected, ngen=ngen)
  
  q<-Q[upper.tri(Q, diag=TRUE)]
  a<-0.5*A[upper.tri(A, diag=TRUE)]
  b<-diag(nrow(Q))[upper.tri(Q, diag=TRUE)]
  
  fit_BC = lm(q~a+b-1)
  
  Ne_BC = 0.5/coef(fit_BC)["b"]
  vA_BC = coef(fit_BC)["a"]
  
  return(list(Ne_BC = Ne_BC, vA_BC = vA_BC))
  
}
