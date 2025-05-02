\name{predict_Ne}
\alias{predict_Ne}
\title{
Predict the Ne of the population
}
\description{
Predict the Ne of the population given environmental variance in fitness only.
}
\usage{
predict_Ne(Ve, n, fitness_model = "Exp")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ve}{Environmental variance in fitness on some scale (log if \code{fitness_model = "Exp"} identity if \code{fitness_model = "I"}  )
}
  \item{n}{Census population size
}
  \item{fitness_model}{Either \code{"Exp"} where fitness is log normal with a meanlog of zero or \code{"I"} where fitness is normal with a mean of one}
}
}
\details{
}
\value{
predicted N_e given environmental variance in fitness only.
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
function (Ve, n, fitness_model = "Exp") 
{
    if (fitness_model == "Exp") {
        Vo <- 4 * exp(Ve) - 2
    }
    if (fitness_model == "I") {
        Vo <- 4 * Ve + 2
    }
    Ne <- 4 * n/(2 + Vo)
    return(Ne)
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
