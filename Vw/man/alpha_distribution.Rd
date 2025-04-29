\name{alpha_distribution}
\alias{alpha_distribution}
\title{
Model fitting for average effects
}
\description{
Function for estimating the parameters of a model that relates average effects to allele frequencies.
}
\usage{
alpha_distribution(alpha, p, tprojp = NULL, logit = FALSE, save_model = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{
A vector of average effects for relative fitness
}
  \item{p}{
A vector of allele frequencies in the parents' generation
}
  \item{tprojp}{
Optional projection matrix (transposed)  - if non-null then analysis carried out on projected allele frequencies
}
  \item{logit}{
Logical: if \code{TRUE} a logit model is used (see details)
}
  \item{save_model}{
Logical: if \code{TRUE} the linear model, conditional on the non-linear parameters (\code{palpha} if \code{logit=FALSE} or \code{palpha} and \code{balpha_0} if \code{logit=FALSE}), is returned. 
}
}
\details{
If \code{logit=FALSE}, the average effects are modelled with mean \code{balpha_0+balpha_1*(p-q)} and variance \code{sigm2alpha*(p*q)^palpha}. If \code{logit=TRUE}, the average effects are modelled with mean \code{balpha_1*logit((sign(p-q) * abs(p-q)^balpha_0 + 1)/2)} and variance \code{sigm2alpha*(p*q)^palpha}.
}
\value{

\item{palpha}{Estimate of \code{palpha}}
\item{balpha_0}{Estimate of \code{balpha_0}}
\item{balpha_1}{Estimate of \code{balpha_1}}
\item{sigma2alpha}{Estimate of \code{sigma2alpha}}
\item{logit}{Logical: if \code{TRUE} a logit model was used}
\item{model}{Linear model object, fitted conditional on the non-linear parameters}
\item{X}{Design matrix of conditional linear model}
\item{S}{Sampling covariance matrix of linear parameters: \code{balpha_0} and \code{balpha_1} if \code{logit=FALSE} or \code{balpha_1} if \code{logit=TRUE}}

}
\references{
}
\author{
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
