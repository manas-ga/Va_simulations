\name{prob_novelA}
\alias{prob_novelA}
\title{
Probability of a novel recombinant sequence.
}
\description{
Calculating the probability that a recombinant sequence produced through gametes pairing up at random according to their frequencies is not present in a random sample of the parental haplotypes.
}
\usage{
prob_novelA(pA, pB, rD, nH)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pA}{Reference allele frequencies at locus A}
  \item{pB}{Reference allele frequencies at locus B}
  \item{rD}{Linkage-disequilibirum measured as a correlation}
  \item{nH}{Number of parental haplotypes}
}

\value{
Probability of a novel recombinant sequence.
}

\references{
Manas Geeta Arun, Aidan Angus-Henry, Darren J. Obbard & Jarrod D. Hadfield (2015) "Estimating the additive genetic variance for relative fitness from changes in allele frequency".
}
\author{
Manas Geeta Arun and Jarrod Hadfield
}