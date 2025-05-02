\name{Vw_model}
\alias{Vw_model}
\title{
Function for estimating the additive genetic variance in fitness
}
\description{
Function for estimating the additive genetic variance in fitness from allele frequency change data.
}
\usage{
Vw_model(c_genome=NULL, nR, pbar0, pbar1, ngen1 = 1, pbar2, ngen2, nind=NULL, proj, LDalpha, pa, palpha, balpha, Vs, method, L, Ltilde, svdL = NULL, Ne = NULL, NE=Ne, tol = sqrt(.Machine$double.eps), save_tprojp = FALSE, verbose = TRUE)
}
\arguments{
  \item{c_genome}{Optional matrix of gametic contributions in the base population with rows being gametes and columns being loci. Elements are one (has the reference allele) or zero (has the alternate allele). Gametes are ordered within individuals (i.e. gametes 1 and 2 are from individual A, gametes 3 and 4 are from individual B)}
  \item{nR}{Matrix of non-recombinant probabilities}
  \item{pbar0}{Vector of allele frequencies in the base population.}
  \item{pbar1}{Matrix of allele frequencies in generation ngen1 with loci in columns and replicate populations in rows.}
  \item{ngen1}{Generation at which allele frequencies are first sampled. Generation 0 is the base population.}
  \item{pbar2}{Matrix of allele frequencies in generation ngen2 with loci in columns and replicate populations in rows.}
  \item{ngen2}{Generation at which allele frequencies are last sampled. Generation 0 is the base population.}
  \item{nind}{Number of individuals in base population. If not given, obtained from c_genome.}
  \item{proj}{Type of projection matrix for allele frequencies. Can be one of "LoM", "BLoM", "L" or "N".}
  \item{LDalpha}{logical, if \code{TRUE} L is considered while modelling distribution of alphas rather than diag(L).}
  \item{pa}{Power to which the singular values are taken when projecting.}
  \item{palpha}{Covariance structure for loci effects are modelled as \code{L^palpha}. If NA, palpha is estimated.}
  \item{balpha}{Mean loci effects are modelled as \code{balpha[1]+(p-q)*balpha[2]}. If NA, balpha is estimated.}
  \item{Vs}{Type of covariance structure for loci effects: "L" or "LoNL".}
  \item{method}{"REML" or "MCMC"}
  \item{L}{Covariance matrix of genotypes}
  \item{Ltilde}{Preserved covariance matrix of genotypes.}
  \item{svdL}{Singular value decomposition of genotype matrix.}
  \item{Ne}{Effective population size: if \code{NULL} set to \code{nind}, if length 1, same Ne assumed throughout, if length 2 different Ne's in the neutral 0:ngen1 (Ne[1]) and selected ngen1:ngen2 (Ne[2]) parts of the experiment, and if a vector of length ngen2 different Ne in each generation.}
  \item{NE}{Effective population size accounting for environmental variation in fitness only: same options as for \code{Ne}.}
  \item{tol}{Value below which singular values are considered to be zero.}
  \item{save_tprojp}{logical, if \code{TRUE} transpose of the projection matrix for allele frequencies is stored}
  \item{verbose}{Should progress be printed to screen?}
}

\value{
  \item{Vw_est}{Estimate of the additive genetic variance for fitness}
  \item{data}{Formatted data for analysis}
  \item{model}{Model object}
  \item{SC}{Covariance structure for locus effects.}
  \item{palpha}{Estimate of palpha}
  \item{balpha}{Estimate of balpha}
  \item{palpha_var}{Sampling variance of palpha}
  \item{balpha_var}{Sampling (co)variance matrix for balpha: S}
  \item{tprojp}{Transpose of the projection matrix for allele frequencies}
  \item{X}{Fixed effect design matrix of model}
  \item{DL}{Square root of the eigenvalues of L}
}
\references{
Manas Geeta Arun, Aidan Angus-Henry, Darren J. Obbard & Jarrod D. Hadfield (2015) "Estimating the additive genetic variance for relative fitness from changes in allele frequency".
}
\author{
Manas Geeta Arun and Jarrod Hadfield
}
