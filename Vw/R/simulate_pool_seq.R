simulate_pool_seq = function(c_genome,              # Matrix of haploid genomes of all individuals with only segregating sites (0 for non-reference allele, 1 for reference allele); rows = genomes - the two genomes of an individual are together, columns = sites 
                    SNPs,                   # positions of segregating sites, in the same order as in c_genome
                    sequence_length,
                    read_length,
                    coverage,
                    OD = 2){                # The ratio of variance:mean (if OD==1 rpois is used to map reads to individuals, if OD>1 rnbinom is used
  
  # Sanity checks
  
  #if(!is.integer(sequence_length)){stop("sequence_length must be an integer")}
  if(ncol(c_genome)>sequence_length){stop("c_genome cannot have more columns (sites) than sequence_length")}
  if(length(SNPs)!=ncol(c_genome)){stop("the length of SNPs must equal the number of columns in c_genome")}
  
  # Calculate expected number of reads to be sampled
  n_reads = round(sequence_length*coverage/read_length, 1)
  n_ind = nrow(c_genome)/2
  n_sites = length(SNPs)
  
  # Permissible range for starting points for reads
  read_range = 1:(sequence_length - read_length + 1)
  
  # Number of reads to be mapped to each individual drawn from a a poisson (if OD = 1) or negative binomial distribution (OD >1)
  
  if(OD == 1){
    ind_reads = rpois(n = n_ind, lambda = n_reads/n_ind)
  }else{
    ind_reads = rnbinom(n = n_ind, mu = n_reads/n_ind, size = (n_reads/n_ind)*(1/(OD-1)))
  }
  
  # Create an empty vector to record how many times each segregating site is hit by a read
  count_total = rep(0, length(SNPs))
  
  # Create an empty vector to record how many times each segregating site is hit by a read and is recorded to have the reference allele
  count_ref = rep(0, length(SNPs))
  
  # Loop over individuals mapping reads
  # For every read note the state of sites within SNPs 
  # Update count_total and count_ref
  
  for(ind in 1:n_ind){
    print(ind)
    # Loop over the total number of reads for this individual (stored within ind_reads)
    
    for (read in 1:ind_reads[ind]){
      # Randomly chose either of the two genomes of this individual (0 or 1)
      genome = sample(c(1,2), size = 1) # This genome will be [2*(n_ind - 1) + genome]th row in c_genome
      
      # Randomly sample a starting position on the genome
      start_pos = sample(read_range, size = 1)
      end_pos = start_pos + read_length - 1
      
      # Identify how many segregating sites are covered by this read
      read_sites = which(SNPs>=start_pos&SNPs<=end_pos)
      
      # Perform subsequent operations only if at least one segregating site is covered by the read
      if(length(read_sites)>0){
        count_total[read_sites] = count_total[read_sites] + 1
        count_ref[read_sites] = count_ref[read_sites] + c_genome[2*(ind - 1) + genome, read_sites] 
      }
      
    }
    
  }
  
  # identify sites that are not covered at all
  
  uncovered_sites = which(count_total == 0)
  
  
  if(length(uncovered_sites>0)){
    warning("Some segregating sites were not covered by any reads")
    count_total[uncovered_sites] = NA
    }
  
  return(list("p" = count_ref/count_total, "coverage" = count_total))
  
  
  
}