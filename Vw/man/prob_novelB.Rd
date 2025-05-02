\name{prob_novelB}
\alias{prob_novelB}
\title{
Probability of a novel recombinant sequence.
}
\description{
Calculating the probability that a recombinant sequence produced through gametes pairing up at random according to their frequencies is not present in a random sample of the parental haplotypes.
}
\usage{
prob_novelB(P, Q, R, S, nH)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{P}{Frequency of the haplotype with reference alleles at both loci}
  \item{Q}{Frequency of the haplotype with reference allele at locus A and alternate allele at locus B}
  \item{R}{Frequency of the haplotype with alternate allele at locus A and reference allele at locus B}
  \item{S}{Frequency of the haplotype with alternate alleles at both loci}
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
