\name{predict_NE}
\alias{predict_NE}
\title{
Predict the NE of the population
}
\description{
Predict the NE of the population given environmental variance in fitness only.
}
\usage{
predict_NE(Ve, n, fitness_model = "Exp")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ve}{Environmental variance in fitness on some scale (log if \code{fitness_model = "Exp"} identity if \code{fitness_model = "I"})}
  \item{n}{Census population size}
  \item{fitness_model}{Either \code{"Exp"} where fitness is log normal with a meanlog of zero or \code{"I"} where fitness is normal with a mean of one}
}

\value{
predicted Ne given environmental variance in fitness only: NE
}

\references{
Manas Geeta Arun, Aidan Angus-Henry, Darren J. Obbard & Jarrod D. Hadfield (2015) "Estimating the additive genetic variance for relative fitness from changes in allele frequency".
}
\author{
Manas Geeta Arun and Jarrod Hadfield
}
