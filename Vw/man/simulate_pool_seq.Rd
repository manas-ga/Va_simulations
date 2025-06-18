\name{simulate_pool_seq}
\alias{simulate_pool_seq}
\title{
Sample allele frequencies from a population by simulating pool-seq.
}
\description{
Read the genotypes of all the genomes within a population, simulate pool-seq on the population, and return sampled allele frequencies and the coverage for each segregating site.
}
\usage{
  simulate_pool_seq(c_genome, pos, sequence_length, read_length, coverage, V_logmean = 0)
}
\arguments{
  \item{c_genome}{Matrix of gametic contributions with rows being gametes and columns being loci. Elements should be one (has the reference allele) or zero (has the alternate allele). Gametes should be ordered
 within individuals (i.e. gametes 1 and 2 are from individual A, gametes 3 and 4 are from individual B).}
  \item{SNPs}{Vector of positions of segregating sites.}
  \item{sequence_length}{Genome size in base-pairs.}
  \item{read_length}{The length of each read in base-pairs.}
  \item{coverage}{How many times on average should each site be covered.}
  \item{V_logmean}{The number of reads mapping to an individual are drawn from a poisson whose means are lognormal; i.e. the logged means are normal with variance = V_logmean (defaults to regular poisson) and mean = log(n_reads/n_ind) - V_logmean/2}
}

\value{
\item{p}{Vector of sampled allele frequencies for the segregating sites.}
\item{coverage}{Vector of the number of reads mapping to the segregating sites.}
}
\references{
Manas Geeta Arun, Aidan Angus-Henry, Darren J. Obbard & Jarrod D. Hadfield (2015) "Estimating the additive genetic variance for relative fitness from changes in allele frequency".
}
\author{
Manas Geeta Arun and Jarrod Hadfield
}

\seealso{
\code{\link{extract_slim_data}}
}
