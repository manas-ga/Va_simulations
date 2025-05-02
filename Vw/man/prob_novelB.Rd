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
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Probability of a novel recombinant sequence.
}
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (P, Q, R, S, nH) 
{
    P * S * (1 - Q)^nH + P * S * (1 - R)^nH + Q * R * (1 - P)^nH + 
        Q * R * (1 - S)^nH
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
