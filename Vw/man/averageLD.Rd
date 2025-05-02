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
  \item{L}{Matrix with elements proportional to the LD between sites in the base population. For the most accurate result Ltilde should be used: see \link{}}
  \item{nR}{Matrix of non-recombination probabilities}
  \item{nrep}{The number of replicate populations}
  \item{selected}{Optional indices for selected sites. If null the LD between all sites is used, if not null the LD between selected and neutral sites is used.}
  \item{ngen}{Number of generations that separate replicates from the base-population. \code{ngen=1} is for a situation where replicates have been initiated from the offspring of the base-population.}
}
\details{
}
\value{
A nrep by nrep matrix with one among the diagonal and the average (over pairs of loci) expected products of linkage-disequilibria (LD) for pairs of replicates o the off-diagonal.
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
function (L, nR, N = Inf, nrep) 
{
    n <- nrow(L)
    s <- (sum(abs(cov2cor(L)) * nR) - n)/(n * (n - 1))
    Sigma <- matrix(s, nrep, nrep)
    diag(Sigma) <- 1 + 1/N
    return(Sigma)
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
