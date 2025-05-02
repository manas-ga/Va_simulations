\name{summaryG}
\alias{summaryG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
\summaryG
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summaryG(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Matrix with two columns with individuals in rows and haplotypes in columns: 1 = reference allele at both loci, 2 = reference allele at locus A and alternate allele at locus B,  3 = alternate allele at locus A and reference allele at locus B and  4  = alternate allele at both loci.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{pA}{Reference allele frequency at locus A}
  \item{pB}{Reference allele frequency at locus A}
  \item{Lgp}{Gametic phase disequilibrium: 2L'}
  \item{Lngp}{Non-gametic phase disequilibrium: 2L''}

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
function (x) 
{
    if (ncol(x) != 2 | any(!x \%in\% c(1:4))) {
        stop("x should be an nx2 matrix with integers 1:4 corresponding to haplotypes AB (1), Ab (2), aB (3) and ab (4)")
    }
    n_hap <- length(x)
    pA <- sum(x == 1 | x == 2)/n_hap
    pB <- sum(x == 1 | x == 3)/n_hap
    Lgp <- sum(x == 1)/n_hap
    Lgp <- Lgp - pA * pB
    Lngp <- sum((x[, 1] == 1 | x[, 1] == 2) & (x[, 2] == 1 | 
        x[, 2] == 3))
    Lngp <- Lngp + sum((x[, 2] == 1 | x[, 2] == 2) & (x[, 1] == 
        1 | x[, 1] == 3))
    Lngp <- Lngp/n_hap
    Lngp <- Lngp - pA * pB
    return(list(pA = pA, pB = pB, Lgp = Lgp, Lngp = Lngp))
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
