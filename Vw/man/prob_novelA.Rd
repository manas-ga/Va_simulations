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
}
  \item{pB}{Reference allele frequencies at locus B}
  \item{rD}{Linkage-disequilibirum measured as a correlation}
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
function (pA, pB, rD, nH) 
{
    pa <- 1 - pA
    pb <- 1 - pB
    D <- rD * sqrt(pA * pa * pB * pb)
    P <- pA * pB + D
    Q <- pA * pb - D
    R <- pa * pB - D
    S <- pa * pb + D
    (P * S * (1 - Q)^nH + P * S * (1 - R)^nH + Q * R * (1 - P)^nH + 
        Q * R * (1 - S)^nH)
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
