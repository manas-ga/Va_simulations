\name{plot_alpha_distribution}
\alias{plot_alpha_distribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot properties of the joint distribution of alphas and allele frequencies
}
\description{
Plot properties of the joint distribution of alphas and allele frequencies
}
\usage{
plot_alpha_distribution(alpha = NULL, p = NULL, tprojp = NULL, parameters, ...)
}
\arguments{
  \item{alpha}{vector of average effects}
  \item{p}{vector of allele frequencies}
  \item{tprojp}{optional projection matrix (transposed) for allele frequencies and average effects}
  \item{parameters}{a list containing the elements \code{logit}, \code{balpha_0}, \code{balpha_1} and \code{palpha}. \code{logit} should be logical indicating whether the model is on the logit scale or not, and the remaining elements should be the numeric values of the parameters: see \code{\link{alpha_distribution}}}
  \item{\dots}{additional arguments passed to the plot function}
}

\references{
Manas Geeta Arun, Aidan Angus-Henry, Darren J. Obbard & Jarrod D. Hadfield (2015) "Estimating the additive genetic variance for relative fitness from changes in allele frequency".
}
\author{
Manas Geeta Arun and Jarrod Hadfield
}
