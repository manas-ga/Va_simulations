
predict_NE<-function(Ve, n, fitness_model = "Exp"){
  
  if(fitness_model == "Exp"){
    # assuming meanlog=0
    Vo<-4*exp(Ve)-2
  }
  if(fitness_model == "I"){
    # assuming mean=1
    Vo<-4*Ve+2
  }
  Ne<-4*n/(2+Vo)
  return(NE)
}