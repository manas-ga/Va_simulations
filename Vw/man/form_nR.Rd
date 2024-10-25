\name{form_nR}
\alias{form_nR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
\Calculate the matrix of non-recombinant fractions
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
form_nR(SNPs, RecombRate, HapLength, AtleastOneRecomb)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SNPs}{
%%     ~~Describe \code{SNPs} here~~
}
  \item{RecombRate}{
%%     ~~Describe \code{RecombRate} here~~
}
  \item{HapLength}{
%%     ~~Describe \code{HapLength} here~~
}
  \item{AtleastOneRecomb}{
%%     ~~Describe \code{AtleastOneRecomb} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
