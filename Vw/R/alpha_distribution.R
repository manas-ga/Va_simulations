
alpha_distribution<-function(alpha, p, tprojp=NULL, logit=FALSE, save_model=FALSE){
  
  pmq<-2*p-1
  pq<-p*(1-p)/2
  
  palpha_fit<-function(palpha, alpha, pmq, pq, tprojp){
    
    if(!is.null(tprojp)){
      alpha<-t(alpha%*%tprojp)
      pmq<-t(pmq%*%tprojp)
      pq<-diag(t(tprojp)%*%(tprojp*pq^palpha))
    }else{
      pq<-pq^palpha
    }
    
    -logLik(lm(alpha~pmq, weights=1/pq))
  }
  
  logit_fit<-function(par, alpha, pmq, pq,  tprojp){
    
    balpha_0<-par[1]
    palpha<-par[2]
    
    pmq<-pmq_trans(pmq, balpha_0)
    
    if(!is.null(tprojp)){
      alpha<-t(alpha%*%tprojp)
      pmq<-t(pmq%*%tprojp)
      pq<-diag(t(tprojp)%*%(tprojp*pq^palpha))
    }else{
      pq<-pq^palpha
    }
    -logLik(lm(alpha~pmq-1, weights=1/pq))
  }
  
  if(logit){
    palpha<-optim(c(1, 1), logit_fit, alpha=alpha, pmq=pmq, pq=pq, tprojp=tprojp)$par
    
    balpha_0<-palpha[1]
    palpha<-palpha[2]
    
    pmq<-pmq_trans(pmq, balpha_0)
    
    if(!is.null(tprojp)){
      alpha<-t(alpha%*%tprojp)
      pmq<-t(pmq%*%tprojp)
      pq<-diag(t(tprojp)%*%(tprojp*pq^palpha))
    }else{
      pq<-pq^palpha
    }
    
    m1<-lm(alpha~pmq-1, weights=1/pq)
  }else{
    palpha<-optim(0, palpha_fit, alpha=alpha, pmq=pmq, pq=pq, tprojp=tprojp, method="Brent", lower=-5, upper=5)$par
    
    if(!is.null(tprojp)){
      alpha<-t(alpha%*%tprojp)
      pmq<-t(pmq%*%tprojp)
      pq<-diag(t(tprojp)%*%(tprojp*pq^palpha))
    }else{
      pq<-pq^palpha
    }
    
    m1<-lm(alpha~pmq, weights=1/pq)
  }  
  
  
  result<-list(palpha=palpha, 
               balpha_0=if(logit){balpha_0}else{coef(m1)["(Intercept)"]}, 
               balpha_1=coef(m1)["pmq"],
               sigma2alpha=summary(m1)$sigma^2,
               logit=logit,
               model=if(save_model){m1}else{NULL},
               X=model.matrix(m1),
               S=summary(m1)$cov.unscaled*summary(m1)$sigma^2)
  
  return(result)
}
