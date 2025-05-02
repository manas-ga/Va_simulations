\name{est_Va_bc}
\alias{est_Va_bc}
\title{
Estimate of the additive genetic variance for fitness using the method of Buffalo and Coop (2019).
}
\description{
Estimate of the additive genetic variance for fitness using the method of Buffalo and Coop (2019).
}
\usage{
est_Va_bc(pbar1, pbar2, L, nR, nrep, selected=NULL, exact=FALSE, ngen=1)
}
\arguments{
  \item{pbar1}{
Matrix of allele frequencies in generation t with loci in columns and replicate populations in rows.
}
  \item{pbar2}{
Matrix of allele frequencies in generation t+1 with loci in columns and replicate populations in rows.
}
  \item{L}{
Matrix with elements proportional to the LD between sites in the population at time t=0. For the most accurate result Ltilde should be used: see \link{}
}
  \item{nrep}{The number of replicate populations}
  \item{selected}{Optional indices for selected sites. If null the LD between all sites is used, if not null the LD between selected and neutral sites is used.}
  \item{ngen}{Number of generations that separate the replicates from the base-population: t. \code{ngen=1} is for a situation where replicates have been initiated from the offspring of the base-population.}
\item{exact}{logical: if \code{TRUE} the original implementation is used, otherwise allele frequencies are projected before taking the covariances in allele frequency chages. 
}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

\item{Ne_BC}{Estimate of the effective population size.}
\item{vA_BC}{Estimate of the additive genetic variation for fitness.}

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
function (pbar1, pbar2, L, nR, nrep) 
{
    Q <- Q(pbar1, pbar2)
    Sigma <- Sigma(L, nR, nrep = nrep)
    q <- Q[upper.tri(Q, diag = TRUE)]
    a <- Sigma[upper.tri(Sigma, diag = TRUE)]
    b <- diag(nrow(Q))[upper.tri(Q, diag = TRUE)]
    coef(lm(q ~ a + b - 1))["a"]
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
