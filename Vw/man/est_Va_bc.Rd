\name{est_Va_bc}
\alias{est_Va_bc}
\title{
Estimate of the additive genetic variance for fitness using the method of Buffalo and Coop (2019).
}
\description{
Estimate of the additive genetic variance for fitness using the method of Buffalo and Coop (2019).
}
\usage{
est_Va_bc(pbar1, pbar2, L, nR, selected=NULL, exact=FALSE, ngen=1)
}
\arguments{
  \item{pbar1}{Matrix of allele frequencies in generation t with loci in columns and replicate populations in rows.}
  \item{pbar2}{Matrix of allele frequencies in generation t+1 with loci in columns and replicate populations in rows.}
  \item{L}{Matrix with elements proportional to the LD between sites in the population at time t=0. For the most accurate result Ltilde should be used}
  \item{nR}{Matrix of non-recombination probabilities}
  \item{selected}{Optional indices for selected sites. If null the LD between all sites is used, if not null the LD between selected and neutral sites is used.}
  \item{ngen}{Number of generations that separate the replicates from the base-population: t. \code{ngen=1} is for a situation where replicates have been initiated from the offspring of the base-population.}
\item{exact}{logical: if \code{TRUE} the original implementation is used, otherwise allele frequencies are projected before taking the covariances in allele frequency chages. 
}

}

\value{

\item{Ne_BC}{Estimate of the effective population size.}
\item{vA_BC}{Estimate of the additive genetic variation for fitness.}

}

\references{
Vince Buffalo & Graham Coop. (2019) "The linked selection signature of rapid adaptation in temporal genomic data." Genetics 213.3 1007-1045.
Manas Geeta Arun, Aidan Angus-Henry, Darren J. Obbard & Jarrod D. Hadfield (2015) "Estimating the additive genetic variance for relative fitness from changes in allele frequency".
}
\author{
Manas Geeta Arun and Jarrod Hadfield
}
