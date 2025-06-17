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
  
  # Create an empty matrix having the same dimensions as c_genome to record how many reads mapped to each segregating site on each genome
  n_mapped_reads  = matrix(0, nrow = nrow(c_genome), ncol = ncol(c_genome))
  
  # Create an empty vector to record how many times each segregating site is hit by a read and is recorded to have the reference allele
  count_ref = rep(0, length(SNPs))
  
  # Loop over individuals mapping reads
  # For every read note the state of sites within SNPs 
  # Update count_total and count_ref
  
  for(ind in 1:n_ind){
    print(ind)
    
    # Decide from a binomial, how many of this individual's reads map to genome 1 and genome 2
    
    n_reads_g1 = rbinom(n = 1, size = ind_reads[ind], p = 0.5)
    n_reads_g2 = ind_reads[ind] - n_reads_g1
    
    # Randomly sample a vector of starting positions on genome 1 and genome 2
    
    start_pos_g1 = sample(read_range, size = n_reads_g1, replace = TRUE)
    start_pos_g2 = sample(read_range, size = n_reads_g2, replace = TRUE)
    
    # Explicitly write down positions covered by each read
    # start_pos_g1 on row1
    # start_pos_g1 + 1 on row 2
    # start_pos_g1 + 2 on row 3
    mapped_pos_g1 = matrix(NA, nrow = read_length, ncol = n_reads_g1)
    mapped_pos_g2 = matrix(NA, nrow = read_length, ncol = n_reads_g2)
    
    for(row in 1:nrow(mapped_pos_g1)){
      mapped_pos_g1[row,] = start_pos_g1 + row - 1
      mapped_pos_g2[row,] = start_pos_g2 + row - 1
    }
    
    mapped_pos_g1 = c(mapped_pos_g1)
    mapped_pos_g2 = c(mapped_pos_g2)

    
    # Retain only those positions which correspond to segregating sites (i.e. are included in SNPs)
    mapped_pos_g1 = mapped_pos_g1[which(mapped_pos_g1%in%SNPs)]
    mapped_pos_g2 = mapped_pos_g2[which(mapped_pos_g2%in%SNPs)]
    
    # Populate n_mapped_reads for the two genomes of the current individual
    
    for(site in 1:n_sites){
      # genome 1
      n_mapped_reads[2*ind - 1, site] = sum(mapped_pos_g1 == SNPs[site])
      # genome 2
      n_mapped_reads[2*ind, site] = sum(mapped_pos_g2 == SNPs[site])
      
    }
    
  }
  
  count_total = colSums(n_mapped_reads)
  
  # identify sites that are not covered at all
  uncovered_sites = which(count_total == 0)
  
  if(length(uncovered_sites>0)){
    warning("Some segregating sites were not covered by any reads")
    count_total[uncovered_sites] = NA
    }
  
  p = colSums(c_genome*n_mapped_reads)/count_total
  
  return(list("p" = p, "coverage" = count_total))
  
  
  
}