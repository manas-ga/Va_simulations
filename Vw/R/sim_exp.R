
sim_exp<-function(c_genome, nind, r, position, nrep, alpha=NULL, Ve=0, fitness_model="Exp"){
  
  if(ncol(c_genome)!=length(position)){stop("position should be the same length as number of columns in c_genome")}
  
  if(!is.null(alpha)){
    if(ncol(c_genome)!=length(alpha)){stop("alpha should be the same length as number of columns in c_genome")}
  }
  
  reordering<-order(position)
  
  c_genome<-c_genome[,reordering]
  position<-position[reordering]
  
  # reorder SNPs according to map position
  
  n0_individuals<-nrow(c_genome)/2
  n_loci<-ncol(c_genome)
  
  pbar<-matrix(NA, nrep, n_loci)
  
  if(!is.null(alpha)){
    alpha<-alpha[reordering]
    genetic_fitness<-c_genome%*%alpha/2
    genetic_fitness<-genetic_fitness[seq(1, 2*n0_individuals, 2)]+genetic_fitness[seq(2, 2*n0_individuals, 2)]
  }else{
    genetic_fitness<-rep(0, n0_individuals)
  }
  
  genetic_fitness<-1+genetic_fitness-mean(genetic_fitness)
  # rescale to have mean of 1 (not necessary under the exp model)  
  
  if(any(genetic_fitness<0)){
    warning(paste(sum(genetic_fitness<0), "genetic fitnesses negative - set to zero"))
    genetic_fitness[which(genetic_fitness<0)]<-0
  }
  
  for(rep in 1:nrep){
    
    fitness<-genetic_fitness
    if(Ve!=0){
      if(fitness_model=="Exp"){
        fitness<-fitness+rnorm(n0_individuals, 0, sqrt(Ve))
      }
      if(fitness_model=="I"){
        m<-min(fitness)
        fitness<-fitness+rgamma(n0_individuals, shape=(m^2)/Ve, scale=Ve/m)-m
        # We wish to generate random environmental variables, X, with E[X]=0,  VAR(X)=Ve and min(X) = -min(genetic_fitness) =-m  so that fitness is non-negative 
        
        # To do so we simulate Y from a gamma and generate X = Y-m to satisfy the constraint that X+G>0 where G is fitness based on genotype.
        
        # Have a=shape and s=scale then E[X] = a*s-m and VAR(X) = a*s^2
        
        # Since E[X]=0 then a*s = m and so VAR(X) = m*s and s=VAR(X)/m=Ve/m
        
        # Since a*s = m then a*Ve/m = m and so a = (m^2)/Ve
      }
    }
    
    if(fitness_model=="Exp"){
      fitness<-exp(fitness)
    }
    
    genomes<-matrix(0, nind, n_loci)
    
    for(i in 1:nind){
      
      for(j in 1:2){  # sample gametes from both parents
        
        parent<-sample(1:n0_individuals, 1, prob=fitness)
        # sample parent j of individual i
        
        gamete<-rep(rbinom(1, prob=0.5, size=1), n_loci)
        # sample grandparental gamete of parent j of individual i
        
        r_event<-reda::simEventData(rho=r, end=max(position))
        # simulate recombination events
        
        if(nrow(r_event)>0){
          r_event<-r_event$time[which(r_event$event==1)]
          for(k in 1:length(r_event)){
            gamete[which(position>r_event[k])]<-abs(1-gamete[which(position>r_event[k])])
          } 
        }
        # switch grandparental gamete when a recombination event happens and so grandparental gamete is a string of 0's and 1's 
        
        gam0<-which(gamete==0)
        # which positions does individual i inherit from grandparental gamete 0 of parent j
        
        if(length(gam0)!=0){     
          genomes[i,gam0]<-genomes[i,gam0]+c_genome[(parent-1)*2+1,gam0]
          # add any reference alleles inherited from grandparental gamete 0 
        }
        if(length(gam0)!=n_loci){
          if(length(gam0)!=0){
            genomes[i,-gam0]<-genomes[i,-gam0]+c_genome[(parent-1)*2+2,-gam0]
          }else{
            genomes[i,]<-genomes[i,]+c_genome[(parent-1)*2+2,]
          }  
          # add any reference alleles inherited from grandparental gamete 1
        }
      }
    }  
    pbar[rep,]<-colMeans(genomes)/2
  }  
  pbar<-pbar[,order(reordering)]
  return(pbar)
}
