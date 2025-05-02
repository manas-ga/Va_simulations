\name{extract_slim_data}
\alias{extract_slim_data}
\title{
Extract useful information from SLiM outputs
}
\description{
Extract useful information from SLiM outputs
}
\usage{
extract_slim_data(Set_ID, sim = 1, ngen2_optional=NULL, unzip = FALSE, slim_output_path, sim_param_path, extract_genomes_path, extract_mut_path, mutations_path, c_matrix_path, n_sample = NULL, randomise = TRUE, delete_temp_files = TRUE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{n_sample}{Number of individuals sampled from the parents' generation (useful if n_ind_exp is large).}
  \item{randomise}{logical: if \code{TRUE}, the reference allele is randomised.}
  \item{delete_temp_files}{logical: if \code{TRUE}, temporary files are deleted.}
  \item{verbose}{logical, if \code{TRUE} progress printed to screen.}
}

\value{
  \item{c_genome}{Matrix of gametic contributions in the base population with rows being gametes and columns being loci. Elements are one (has the reference allele) or zero (has the alternate allele). Gametes are ordered within individuals (i.e. gametes 1 and 2 are from individual A, gametes 3 and 4 are from individual B)}
  \item{list_alpha}{vector of average effects}
  \item{SNPs}{Vector of loci positions.}
  \item{ngen1}{Generation at which allele frequencies are first sampled. Generation 0 is the base population.}
  \item{ngen2}{Generation at which allele frequencies are last sampled. Generation 0 is the base population.}
  \item{pbar0}{Vector of allele frequencies in the base population.}
  \item{pbar1}{Matrix of allele frequencies in generation ngen1 with loci in columns and replicate populations in rows.}
  \item{pbar2}{Matrix of allele frequencies in generation ngen2 with loci in columns and replicate populations in rows.}
  \item{sim_params}{List containing information about the parameters used in the simulation.}
}

\references{
Manas Geeta Arun, Aidan Angus-Henry, Darren J. Obbard & Jarrod D. Hadfield (2015) "Estimating the additive genetic variance for relative fitness from changes in allele frequency".
}
\author{
Manas Geeta Arun and Jarrod Hadfield
}
