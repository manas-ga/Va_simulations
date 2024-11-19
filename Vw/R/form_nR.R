form_nR<-function(SNPs,                   # SNP positions
                  RecombRate,             # Recombination rate between adjacent sites
                  HapLength,              # Genome size
                  AtleastOneRecomb        # Should there be at least one cross-over per meiosis?
){
  
  nsnps<-length(SNPs)
  
  message("Computing the distance matrix...")
  
  Dist<-as.matrix(dist(SNPs, diag=TRUE, upper=TRUE))
  # conditional on the total number of recombination events (tn) along a chromosome, the number that fall in the interval is binomial with probability equal to the relative length of the interval to the whole and number of trials equal to tn. The non-recombination fraction is then the probability that an even number of events (including zero) fall in the interval. 
  
  message("Computing nR...")
  
  if(AtleastOneRecomb){
    
    maxn<-actuar::qztpois(1-1e-5, RecombRate*HapLength)
    # maximum number of total recombination events likely to be seen
    
    nR<-matrix(0, nsnps, nsnps)
    
    for(nt in 0:maxn){
      # could iterate nt from 1 since Pr(0 recombination events in total)=0.
      for(ne in seq(0,maxn,2)){
        nR<-nR+dbinom(ne, prob=Dist/HapLength, size=nt)*actuar::dztpois(nt, RecombRate*HapLength)
        # Pr(even number of crossovers falling between snps given nt)*Pr(nt recombination events in total) then sum over nt
      }
    }
    nR<-nR/actuar::pztpois(maxn, RecombRate*HapLength)
    # rescale so diag(R)=1
    
  }else{
    
    nR<-1 - 0.5*(1 - exp(-2*Dist*RecombRate))
  }
  
  return(nR)
  
}
