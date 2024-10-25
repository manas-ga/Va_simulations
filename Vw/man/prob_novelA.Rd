\name{prob_novelA}
\alias{prob_novelA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
prob_novelA(pA, pB, rD, nH)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pA}{
%%     ~~Describe \code{pA} here~~
}
  \item{pB}{
%%     ~~Describe \code{pB} here~~
}
  \item{rD}{
%%     ~~Describe \code{rD} here~~
}
  \item{nH}{
%%     ~~Describe \code{nH} here~~
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
