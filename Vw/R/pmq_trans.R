pmq_trans<-function(pmq, balpha_0){
  x<-sign(pmq)*abs(pmq)^balpha_0
  logit((x+1)/2)
}
