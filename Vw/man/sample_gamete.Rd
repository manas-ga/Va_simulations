\name{sample_gamete}
\alias{sample_gamete}
\title{
Samples gametic 2-locus haplotypes from an individual
}
\description{
Samples gametic 2-locus haplotypes from an individual
}
\usage{
sample_gamete(xi, r, n=1)
}
\arguments{
  \item{xi}{Vector of length two indicating the individual's two haplotypes: 1 = reference allele at both loci, 2 = reference allele at locus A and alternate allele at locus B,  3 = alternate allele at locus A and reference allele at locus B and  4  = alternate allele at both loci.}
  \item{r}{Recombination rate between the two loci}
  \item{n}{Number of gametes to sample}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Gametic 2-locus haplotype: 1 = reference allele at both loci, 2 = reference allele at locus A and alternate allele at locus B,  3 = alternate allele at locus A and reference allele at locus B and  4  = alternate allele at both loci.}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
table(sample_gamete(c(1,4), r=1/4, n=10000))/10000
# AA with BB and r=1/4 expect  AA (0.375) BB (0.375), Ab (0.125) and aB (0.125)
}
