\name{Q}
\alias{Q}
\title{
Allele frequency change covariances 
}
\description{
Function to obtain (scaled) allele frequency change covariances (Q) from Buffalo and Coop (2019)
}
\usage{
Q(pbar1, pbar2, exact=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{pbar1}{Matrix of allele frequencies in generation t with loci in columns and replicate populations in rows.}
\item{pbar2}{Matrix of allele frequencies in generation t+1 with loci in columns and replicate populations in rows.}
\item{exact}{logical: if \code{TRUE} the original implementation is used, otherwise allele frequencies are projected before taking the covariances in allele frequency changes.} 

}

\value{
Matrix of (scaled) allele frequency change covariances. 
}
\references{
Vince Buffalo & Graham Coop. (2019) "The linked selection signature of rapid adaptation in temporal genomic data." Genetics 213.3 1007-1045.
Manas Geeta Arun, Aidan Angus-Henry, Darren J. Obbard & Jarrod D. Hadfield (2015) "Estimating the additive genetic variance for relative fitness from changes in allele frequency".
}
\author{
Manas Geeta Arun and Jarrod Hadfield
}
