\name{averageLD}
\alias{averageLD}
\title{
Average persistent associations between sites across replicate popalations.
}
\description{
Calculates the average (over pairs of loci) expected products of linkage-disequilibria (LD) for pairs of replicates.  
}
\usage{
averageLD(L, nR, nrep, selected=NULL, ngen=1)
}
\arguments{
  \item{L}{Matrix with elements proportional to the LD between sites in the base population. For the most accurate result Ltilde should be used}
  \item{nR}{Matrix of non-recombination probabilities}
  \item{nrep}{The number of replicate populations}
  \item{selected}{Optional indices for selected sites. If null the LD between all sites is used, if not null the LD between selected and neutral sites is used.}
  \item{ngen}{Number of generations that separate replicates from the base-population. \code{ngen=1} is for a situation where replicates have been initiated from the offspring of the base-population.}
}

\value{
A nrep by nrep matrix with one among the diagonal and the average (over pairs of loci) expected products of linkage-disequilibria (LD) for pairs of replicates o the off-diagonal.
}

\references{
Vince Buffalo & Graham Coop. (2019) "The linked selection signature of rapid adaptation in temporal genomic data." Genetics 213.3 1007-1045.
Manas Geeta Arun, Aidan Angus-Henry, Darren J. Obbard & Jarrod D. Hadfield (2015) "Estimating the additive genetic variance for relative fitness from changes in allele frequency".
}
\author{
Manas Geeta Arun and Jarrod Hadfield
}
