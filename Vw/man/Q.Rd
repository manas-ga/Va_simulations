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
\item{exact}{logical: if \code{TRUE} the original implementation is used, otherwise allele frequencies are projected before taking the covariances in allele frequency chages. 

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Matrix of (scaled) allele frequency change covariances. 
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
function (pbar1, pbar2) 
{
    if (any(pbar2 == 0) | any(pbar2 == 1)) {
        pbar2[which(pbar2 == 0 | pbar2 == 1)] <- NA
    }
    if (any(pbar1 == 0) | any(pbar1 == 1)) {
        pbar1[which(pbar1 == 0 | pbar1 == 1)] <- NA
    }
    Q <- cov(t(pbar2 - pbar1), use = "pairwise.complete.obs")
    d <- rowMeans(pbar1 * (1 - pbar1), na.rm = TRUE)
    Q <- Q/outer(d, d)
    return(Q)
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
