
plot_alpha_distribution<-function(alpha=NULL, p=NULL, tprojp=NULL, parameters, ...){
  
  pmq<-2*p-1
  pq<-p*(1-p)/2
  
  if(parameters$logit){
    pred<-pmq_trans(pmq, parameters$balpha_0)
    intercept<-0
  }else{
    pred<-pmq
    intercept<-parameters$balpha_0
  }
  slope<-parameters$balpha_1
  
  palpha<-parameters$palpha
  var<-pq^palpha
  
  if(is.null(tprojp)){
    
    pmq.x<-seq(min(pmq),max(pmq), length=1000)
    pq.x = (pmq.x+1)*(1-pmq.x)/8
    sd.x<-parameters$sigma2alpha*sqrt(pq.x^palpha)
    
    plot(alpha~pmq, ylab="alpha", xlab="p-q", ...)
    
    if(parameters$logit){
      pred.x<-pmq_trans(pmq.x, parameters$balpha_0)
    }else{
      pred.x<-pmq.x
    } 
    
    m.x<-intercept+slope*pred.x
    u.x<-qnorm(0.975, intercept+slope*pred.x, sd.x)
    l.x<-qnorm(0.025, intercept+slope*pred.x, sd.x)
    
    lines(m.x~pmq.x, col="red")
    lines(u.x~pmq.x, col="red", lty=2)
    lines(l.x~pmq.x, col="red", lty=2)
    
  }else{  
    
    alpha<-t(alpha%*%tprojp)
    pred<-t(pred%*%tprojp)
    var<-diag(t(tprojp)%*%(tprojp*pq^palpha))
    
    plot(alpha~pred, ylab="projected alpha", xlab="projected predictor", ...)
    
    pred.x<-seq(min(pred),max(pred), length=1000)
    
    m.x<-intercept+slope*pred.x
    lines(m.x~pred.x, col="red")
    
  }
}  