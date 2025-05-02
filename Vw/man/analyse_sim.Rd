\name{analyse_sim}
\alias{analyse_sim}
\title{
Analyse a SLiM simulation. 
}
\description{
Performs a set of analyses on a SLiM simulation including estimating Va for fitness. 
}
\usage{
analyse_sim(Set_ID, sim = 1, unzip = FALSE, slim_output_path, sim_param_path, extract_genomes_path, extract_mut_path, mutations_path, c_matrix_path, output_path, n_sample = NULL, randomise = TRUE, delete_temp_files = TRUE, proj = "BLoM", LDalpha = FALSE, pa = 1, Vs = "LoNL", method = "REML", palpha = NA, balpha = c(NA, NA), AtleastOneRecomb = FALSE, all.gp = FALSE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Set_ID}{Unique ID of the set of simulations that are controlled by a single R script.}
  \item{sim}{integer: the simulation within the set to be analysed.}
  \item{unzip}{logical: if \code{TRUE}, the SLiM output is unzipped, read, and then rezipped}
  \item{slim_output_path}{File path to directory where the SLiM outputs (for parents and experimental replicates) are stored (as .txt files).}
  \item{sim_param_path}{File path to directory where the .csv file containing simulation parameters is stored.}
  \item{extract_genomes_path}{File path to the python script that extracts genomes and mutations from SLim outputs.}
  \item{extract_mut_path}{File path to the python script that extracts mutations from SLim outputs.}
  \item{mutations_path}{File path to directory where extracted mutations are to be stored (temp files).}
  \item{c_matrix_path}{File path to directory where extracted genomes are to be stored (temp files).}
  \item{output_path}{File path to directory where the final data file is to be stored.}
  \item{n_sample}{Number of individuals sampled from the parents' generation (useful if n_ind_exp is large).}
  \item{randomise}{logical: if \code{TRUE}, the reference allele is randomised.}
  \item{delete_temp_files}{logical: if \code{TRUE}, temporary files are deleted.}
  \item{proj}{projection type for allele frequencies: "LoM", "BLoM", "L" or "N".}
  \item{LDalpha}{logical, if \code{TRUE} L is considered while modelling distribution of alphas rather than diag(L).}
  \item{pa}{Power to which the singular values are taken when projecting.}
  \item{Vs}{Type of covariance structure for loci effects: "L" or "LoNL".}
  \item{method}{Model fitted using "REML" or "MCMC".}
  \item{palpha}{Covariance structure for loci effects are modelled as \code{L^palpha}. If NA, palpha is estimated using optim().}
  \item{balpha}{Mean loci effects are modelled as \code{balpha[1]+(p-q)*balpha[2]}. If NA, balpha is estimated.}
  \item{AtleastOneRecomb}{When calculating nR should it be assumed there is a minimum  of one cross-over per meiosis?}
  \item{all.gp}{logical, if \code{TRUE} non-gametic phase disequilibrium is assumed to be zero when calculating Ltilde.} 
  \item{verbose}{logical, if \code{TRUE} progress printed to screen.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A matrix of 
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

}
