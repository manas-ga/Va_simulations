\name{pmq_trans}
\alias{pmq_trans}
\title{
Transforms p-q for use in the logit model.
}
\description{
Transforms p-q for use in the logit model for the joint distribution of alphas and allele frequencies
}
\usage{
pmq_trans(pmq, balpha_0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pmq}{vector of p-q}
  \item{balpha_0}{power parameter: only values greater than one make biological sense.}
}

\value{
\code{logit((sign(p-q) * abs(p-q)^balpha_0 + 1)/2)}
}

\references{
Manas Geeta Arun, Aidan Angus-Henry, Darren J. Obbard & Jarrod D. Hadfield (2015) "Estimating the additive genetic variance for relative fitness from changes in allele frequency".
}
\author{
Manas Geeta Arun and Jarrod Hadfield
}

\examples{
pmq<-seq(-1, 1, length=100)
balpha_0<-seq(1, 10, length=10)
plot(pmq_trans(pmq, balpha_0=balpha_0[1]), type="l", xlab="p-q", ylab = "pmq_trans(p-q)")
for(i in 2:10){
  lines(pmq_trans(pmq, balpha_0=balpha_0[i]))
}
}

