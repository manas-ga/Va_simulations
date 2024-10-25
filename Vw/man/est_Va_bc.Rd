\name{est_Va_bc}
\alias{est_Va_bc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
\est_Va_bc
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
est_Va_bc(pbar1, pbar2, L, nR, nrep)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pbar1}{
%%     ~~Describe \code{pbar1} here~~
}
  \item{pbar2}{
%%     ~~Describe \code{pbar2} here~~
}
  \item{L}{
%%     ~~Describe \code{L} here~~
}
  \item{nR}{
%%     ~~Describe \code{nR} here~~
}
  \item{nrep}{
%%     ~~Describe \code{nrep} here~~
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
