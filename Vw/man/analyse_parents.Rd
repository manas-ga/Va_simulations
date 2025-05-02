\name{analyse_parents}
\alias{analyse_parents}
\title{
Obtain useful information about the parents' generation
}
\description{
Obtain useful information about the parents' generation
}
\usage{
analyse_parents(c_genome, list_alpha, compute_svdL = FALSE, LDalpha = FALSE, SNPs, RecombRate, HapLength, tol = sqrt(.Machine$double.eps), AtleastOneRecomb = FALSE, calc_nR = TRUE, verbose = TRUE)
}
\arguments{
  \item{c_genome}{Matrix of gametic contributions with rows being gametes and columns being loci. Elements should be one (has the reference allele) or zero (has the alternate allele). Gametes should be ordered within individuals (i.e. gametes 1 and 2 are from individual A, gametes 3 and 4 are from individual B)}
  \item{list_alpha}{vector of average effects}
  \item{compute_svdL}{Should the singular value decomposition of c be computed? c is a matrix of genotypes formed from c_gemome with elements 0 (alternate homozygote), 1/2 (heterozygote) or 1 (reference homozygote). If \code{LDalpha=TRUE} singular value decomposition is performed irrespective of \code{compute_svdL}.)}
  \item{LDalpha}{If \code{TRUE}, the covariance of the average effects is proportional to \code{L^palpha} where \code{L=cov(c)}. If \code{FALSE}, \code{L} is diagonalised before performing \code{L^palpha}.}
  \item{SNPs}{Vector of loci positions.}
  \item{RecombRate}{Recombination rate between adjacent sites.}
  \item{HapLength}{Genome size in base-pairs.}
  \item{tol}{Value below which singular values are considered to be zero.}
  \item{AtleastOneRecomb}{When calculating nR should it be assumed there is a minimum  of one cross-over per meiosis?}
  \item{calc_nR}{Should \code{nR}, a matrix of non-recombinant probabilities between loci, be calculated?}
  \item{verbose}{Should progress be printed to screen?}
}

\value{
\item{L}{Covariance matrix of genotypes}
\item{Ltilde}{Preserved covariance matrix of genotypes.}
\item{nR}{Matrix of probabilities of no recombination.}
\item{svdL}{Singular value decomposition of genotype matrix.}
\item{seg_sites}{Number of segregating sites.}
\item{seg_sites_neu}{Number of segregating neutral sites.}
\item{seg_sites_ben}{Number of segregating beneficial sites.}
\item{seg_sites_del}{Number of segregating deleterious sites.}
\item{mean_diversity}{Genome-wide genetic diversity: pq/2, rather than the more usual 2pq.}
\item{theta}{Watterson's theta.}
\item{va_true}{True additive genic variance.}
\item{vA_true}{True additive genetic variance.}
\item{vA_alpha_emp}{Estimate of the additive genetic variance using a model of the average effects on genetic diversity: see \code{\link{alpha_distribution}}}
\item{parameters}{Parameter estimates from a model of the average effects on genetic diversity: see \code{\link{alpha_distribution}}}
}
\references{
Manas Geeta Arun, Aidan Angus-Henry, Darren J. Obbard & Jarrod D. Hadfield (2015) "Estimating the additive genetic variance for relative fitness from changes in allele frequency".
}
\author{
Manas Geeta Arun and Jarrod Hadfield
}

\seealso{
\code{\link{alpha_distribution}}
}
