\name{form_nR}
\alias{form_nR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the matrix of non-recombinant probabilities
}
\description{
Calculates the matrix of non-recombinant probabilities
}
\usage{
form_nR(SNPs, RecombRate, HapLength, AtleastOneRecomb)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SNPs}{Vector of loci positions.}
  \item{RecombRate}{Recombination rate between adjacent sites.}
  \item{HapLength}{Genome size in base-pairs.}
  \item{AtleastOneRecomb}{When calculating nR should it be assumed there is a minimum  of one cross-over per meiosis?}
}

\value{
matrix of non-recombinant probabilities
}

\references{
Manas Geeta Arun, Aidan Angus-Henry, Darren J. Obbard & Jarrod D. Hadfield (2015) "Estimating the additive genetic variance for relative fitness from changes in allele frequency".
}
\author{
Manas Geeta Arun and Jarrod Hadfield
}
