\name{analyse_sim}
\alias{analyse_sim}
\title{
Analyse a SLiM simulation. 
}
\description{
Performs a set of analyses on a SLiM simulation including estimating Va for fitness. 
}
\usage{
analyse_sim(Set_ID, sim = 1, ngen2_optional=NULL, unzip = FALSE, slim_output_path, 
  sim_param_path, extract_genomes_path, extract_mut_path, mutations_path, 
  c_matrix_path, output_path, n_sample = NULL, randomise = TRUE, 
  delete_temp_files = TRUE, pool_seq = FALSE, read_length = NULL, coverage = NULL, V_
logmean = NULL, proj = "BLoM", LDalpha = FALSE, pa = 1, Vs = "LoNL", 
  method = "REML", palpha = NA, balpha = c(NA, NA), AtleastOneRecomb = FALSE, 
  NE=NULL, Ne=NE, predict_NE = FALSE, all.gp = FALSE, verbose = TRUE)
}
\arguments{
  \item{Set_ID}{Unique ID of the set of simulations that are controlled by a single R script.}
  \item{sim}{integer: the simulation within the set to be analysed.}
  \item{ngen2_optional}{Optional generation at which final allele frequencies are obtained. Otherwise the final generation is used.}
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
  \item{pool_seq}{logical, if \code{TRUE} allele frequencies in the experiment are sampled by simulating pool_seq using s
imulate_pool_seq()}
  \item{read_length}{integer, the length of reads to be mapped in simulate_pool_seq()}
  \item{coverage}{integer, the expected number of reads mapping to every site in simulate_pool_seq()}
  \item{V_logmean}{the variance of log means of poisson distributions used to model the number of reads mapping to indiv
iduals in simulate_pool_seq(). The log means are modeled to be normal with variance V_logmean}
  \item{proj}{projection type for allele frequencies: "LoM", "BLoM", "L" or "N".}
  \item{LDalpha}{logical, if \code{TRUE} L is considered while modelling distribution of alphas rather than diag(L).}
  \item{pa}{Power to which the singular values are taken when projecting.}
  \item{Vs}{Type of covariance structure for loci effects: "L" or "LoNL".}
  \item{method}{Model fitted using "REML" or "MCMC".}
  \item{palpha}{Covariance structure for loci effects are modelled as \code{L^palpha}. If NA, palpha is estimated using optim().}
  \item{balpha}{Mean loci effects are modelled as \code{balpha[1]+(p-q)*balpha[2]}. If NA, balpha is estimated.}
  \item{AtleastOneRecomb}{When calculating nR should it be assumed there is a minimum  of one cross-over per meiosis?}
  \item{NE}{Effective population size accounting for environmental variation in fitness only: if \code{NULL} set to \code{nind}, if length 1, same NE assumed throughout, if length 2 different NE's in the neutral 0:ngen1 (NE[1]) and selected ngen1:ngen2 (NE[2]) parts of the experiment, and if a vector of length ngen2 different NE in each generation.}
  \item{Ne}{Effective population size: same options as for \code{NE}.}
  \item{predict_NE}{logical: if \code{TRUE} Ne and NE computed using \code{predict_NE(c(0,Ve_w_expt), n=nind_expt)}}
  \item{all.gp}{logical, if \code{TRUE} non-gametic phase disequilibrium is assumed to be zero when calculating Ltilde.} 
  \item{verbose}{logical, if \code{TRUE} progress printed to screen.}
}
\value{
A matrix of 
}
\references{
Manas Geeta Arun, Aidan Angus-Henry, Darren J. Obbard & Jarrod D. Hadfield (2015) "Estimating the additive genetic variance for relative fitness from changes in allele frequency".
}
\author{
Manas Geeta Arun and Jarrod Hadfield
}
\seealso{
\code{\link{alpha_distribution}}
\code{\link{analyse_parents}}
\code{\link{Vw_model}}
}

