\name{summaryG}
\alias{summaryG}
\title{
Summarises genetic properties of a population from haplotypes
}
\description{
Summarises genetic properties of a population from haplotypes
}
\usage{
summaryG(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Matrix with two columns with individuals in rows and haplotypes in columns: 1 = reference allele at both loci, 2 = reference allele at locus A and alternate allele at locus B,  3 = alternate allele at locus A and reference allele at locus B and  4  = alternate allele at both loci.
}
}

\value{
  \item{pA}{Reference allele frequency at locus A}
  \item{pB}{Reference allele frequency at locus A}
  \item{Lgp}{Gametic phase disequilibrium: 2L'}
  \item{Lngp}{Non-gametic phase disequilibrium: 2L''}

}
\references{
Manas Geeta Arun, Aidan Angus-Henry, Darren J. Obbard & Jarrod D. Hadfield (2015) "Estimating the additive genetic variance for relative fitness from changes in allele frequency".
}
\author{
Manas Geeta Arun and Jarrod Hadfield
}

