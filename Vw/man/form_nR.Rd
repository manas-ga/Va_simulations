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
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
matrix of non-recombinant probabilities
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
function (SNPs, RecombRate, HapLength, AtleastOneRecomb) 
{
    nsnps <- length(SNPs)
    message("Computing the distance matrix...")
    Dist <- as.matrix(dist(SNPs, diag = TRUE, upper = TRUE))
    message("Computing nR...")
    if (AtleastOneRecomb) {
        maxn <- actuar::qztpois(1 - 1e-05, RecombRate * HapLength)
        nR <- matrix(0, nsnps, nsnps)
        for (nt in 0:maxn) {
            for (ne in seq(0, maxn, 2)) {
                nR <- nR + dbinom(ne, prob = Dist/HapLength, 
                  size = nt) * actuar::dztpois(nt, RecombRate * 
                  HapLength)
            }
        }
        nR <- nR/actuar::pztpois(maxn, RecombRate * HapLength)
    }
    else {
        nR <- 1 - 0.5 * (1 - exp(-2 * Dist * RecombRate))
    }
    return(nR)
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
