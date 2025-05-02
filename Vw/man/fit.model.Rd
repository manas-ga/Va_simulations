\name{fit.model}
\alias{fit.model}
\title{
Fit a linear mixed model to allele frequency changes conditional on palpha.
}
\description{
Fit a linear mixed model to allele frequency changes conditional on palpha.
}
\usage{
fit.model(palpha, balpha, LDalpha, nsnps, UL, DL, L, ngen2, ngen1, tprojp, pbar0, pbar1, pbar2, nrep, Selec, LLonly = FALSE, method = "REML", verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{palpha}{Covariance structure for loci effects are modelled as \code{L^palpha}.}
  \item{balpha}{Mean loci effects are modelled as \code{balpha[1]+(p-q)*balpha[2]}. If NA, balpha is estimated.}
  \item{LDalpha}{logical, if \code{TRUE} L is considered while modelling distribution of alphas rather than diag(L).}
  \item{nsnps}{Number of loci}
  \item{UL}{Eigenvectors of L}
  \item{DL}{Eigenvectors of L square-rooted}
  \item{L}{Matrix with elements proportional to the LD between sites in the base population. For the most accurate result Ltilde should be used}
  \item{ngen2}{Generation at which allele frequencies are last sampled. Generation 0 is the base population.}
  \item{ngen1}{Generation at which allele frequencies are first sampled. Generation 0 is the base population.}
  \item{tprojp}{Transpose of the projection matrix for allele frequencies.}
  \item{pbar0}{Vector of allele frequencies in the base population.}
  \item{pbar1}{Matrix of allele frequencies in generation ngen1 with loci in columns and replicate populations in rows.}
  \item{pbar2}{Matrix of allele frequencies in generation ngen2 with loci in columns and replicate populations in rows.}
  \item{nrep}{The number of replicate populations}
  \item{Selec}{Matrix, when post-multiplied by the average effects, gives the expected allele frequency change: Lmathcal}
  \item{LLonly}{logical: if \code{TRUE} only the (restricted) likelihood is returned}
  \item{method}{"REML" or "MCMC"}
  \item{verbose}{logical, if \code{TRUE} progress printed to screen.}
}

\value{
\item{data}{Formatted data for analysis.}
\item{model}{Model object.}
\item{SC}{Covariance structure for locus effects.}
\item{palpha}{Value of palpha used}
%% ...
}

\references{
Manas Geeta Arun, Aidan Angus-Henry, Darren J. Obbard & Jarrod D. Hadfield (2015) "Estimating the additive genetic variance for relative fitness from changes in allele frequency".
}
\author{
Manas Geeta Arun and Jarrod Hadfield
}

