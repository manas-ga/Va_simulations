\name{extract_slim_data}
\alias{extract_slim_data}
\title{
Extract useful information from SLiM outputs
}
\description{
Extract useful information from SLiM outputs
}
\usage{
extract_slim_data(Set_ID, sim = 1, unzip = FALSE, slim_output_path, sim_param_path, extract_genomes_path, extract_mut_path, mutations_path, c_matrix_path, n_sample = NULL, randomise = TRUE, delete_temp_files = TRUE, verbose = TRUE)
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
  \item{n_sample}{Number of individuals sampled from the parents' generation (useful if n_ind_exp is large).}
  \item{randomise}{logical: if \code{TRUE}, the reference allele is randomised.}
  \item{delete_temp_files}{logical: if \code{TRUE}, temporary files are deleted.}
  \item{verbose}{logical, if \code{TRUE} progress printed to screen.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Set_ID, sim = 1, unzip = FALSE, slim_output_path, sim_param_path, 
    extract_genomes_path, extract_mut_path, mutations_path, c_matrix_path, 
    n_sample = NULL, randomise = TRUE, delete_temp_files = TRUE, 
    verbose = TRUE) 
{
    sim_params = read.csv(paste(sim_param_path, "/", Set_ID, 
        "_Data.csv", sep = ""), header = T)
    if ("sim" \%in\% colnames(sim_params)) {
        sim_params = sim_params[sim_params$sim == sim, ]
    }
    if (nrow(sim_params) == 0) {
        stop("No data found!")
    }
    n_ind_exp = sim_params$n_ind_exp
    n_cages = sim_params$n_cages
    if ("ngen1" \%in\% colnames(sim_params)) {
        ngen1 = sim_params$ngen1
    }
    else {
        ngen1 = 1
    }
    if ("ngen2" \%in\% colnames(sim_params)) {
        ngen2 = sim_params$ngen2
    }
    else {
        ngen2 = sim_params$ngen_expt + 1
    }
    ngen_expt = ngen2 - ngen1
    end_gen = sim_params$end_gen
    if (is.null(n_sample)) {
        n_sample = n_ind_exp
    }
    if (n_sample > n_ind_exp) {
        stop("n_sample cannot be greater than n_ind_exp")
    }
    if (verbose) {
        message("Reading the state of the population in the parent's generation...")
        message("Extracting mutations and genomes from the parents' generation...")
    }
    if (unzip) {
        if (verbose) 
            message("Unzipping the SLiM output file...")
        system(paste("gunzip", paste(slim_output_path, "/", Set_ID, 
            "_sim", sim, "_output_parents.txt.gz", sep = "")))
    }
    system(paste("python", extract_genomes_path, paste(slim_output_path, 
        "/", Set_ID, "_sim", sim, "_output_parents.txt", sep = ""), 
        paste(mutations_path, "/", Set_ID, "_sim", sim, "_mutations_parents.txt", 
            sep = ""), paste(c_matrix_path, "/", Set_ID, "_sim", 
            sim, "_c_matrix_parents.csv", sep = ""), n_sample))
    if (unzip) {
        if (verbose) 
            message("Re-zipping the unzipped SLiM output file...")
        system(paste("gzip", paste(slim_output_path, "/", Set_ID, 
            "_sim", sim, "_output_parents.txt", sep = "")))
    }
    c_genome = read.csv(paste(c_matrix_path, "/", Set_ID, "_sim", 
        sim, "_c_matrix_parents.csv", sep = ""), header = F)
    if (delete_temp_files) {
        system(paste("rm", paste(c_matrix_path, "/", Set_ID, 
            "_sim", sim, "_c_matrix_parents.csv", sep = ""), 
            sep = " "))
    }
    c_genome = as.matrix(c_genome)
    n0_individuals = nrow(c_genome)/2
    retained_loci = which(colSums(c_genome) != 0)
    c_genome = c_genome[, retained_loci]
    n_sites = ncol(c_genome)
    pbar0 = colMeans(c_genome)
    mutations_0 = read.table(paste(mutations_path, "/", Set_ID, 
        "_sim", sim, "_mutations_parents.txt", sep = ""), header = T, 
        sep = " ")
    if (delete_temp_files) {
        system(paste("rm ", mutations_path, "/", Set_ID, "_sim", 
            sim, "_mutations_parents.txt", sep = ""))
    }
    mutations_0 = mutations_0[order(mutations_0$Temp_ID), ]
    mutations_0 = mutations_0[retained_loci, ]
    list_alpha = 2 * (mutations_0$s)
    SNPs = mutations_0$Position
    d_proj = c()
    d_raw = c()
    P_matrix = c()
    for (cage in (1:n_cages)) {
        if (verbose) {
            message(paste("Extracting mutations and storing allele frequencies in cage ", 
                cage, " of simulation ", sim, "...", sep = ""))
        }
        for (gen in (end_gen + 1):(end_gen + ngen2)) {
            if (unzip) {
                if (verbose) {
                  message("Unzipping the SLiM output file...")
                }
                system(paste("gunzip", paste(slim_output_path, 
                  "/", Set_ID, "_sim", sim, "_cage", cage, "_output_experiment_", 
                  as.integer(gen), ".txt", sep = "")))
            }
            system(paste("python", extract_mut_path, paste(slim_output_path, 
                "/", Set_ID, "_sim", sim, "_cage", cage, "_output_experiment_", 
                as.integer(gen), ".txt", sep = ""), paste(mutations_path, 
                "/", Set_ID, "_sim", sim, "_cage", cage, "_mutations_", 
                as.integer(gen), ".txt", sep = "")))
            if (unzip) {
                if (verbose) {
                  message("Re-zipping the unzipped SLiM output file...")
                }
                system(paste("gzip", paste(slim_output_path, 
                  "/", Set_ID, "_sim", sim, "_cage", cage, "_output_experiment_", 
                  as.integer(gen), ".txt", sep = "")))
            }
        }
        P = c()
        P = cbind(P, mutations_0$Number/(2 * n_ind_exp))
        for (gen in (end_gen + 1):(end_gen + ngen2)) {
            mut = read.csv(paste(mutations_path, "/", Set_ID, 
                "_sim", sim, "_cage", cage, "_mutations_", as.integer(gen), 
                ".txt", sep = ""), sep = " ")
            if (delete_temp_files) {
                system(paste("rm ", mutations_path, "/", Set_ID, 
                  "_sim", sim, "_cage", cage, "_mutations_", 
                  as.integer(gen), ".txt", sep = ""))
            }
            mut = mut[order(mut$Temp_ID), ]
            freq = c()
            for (mutation in mutations_0$Permanent_ID) {
                if (mutation \%in\% mut$Permanent_ID) {
                  freq = c(freq, (mut[which(mut$Permanent_ID == 
                    mutation), ]$Number)/(2 * n_ind_exp))
                }
                else {
                  freq = c(freq, round(P[which(mutations_0$Permanent_ID == 
                    mutation), gen - end_gen]))
                }
            }
            P = cbind(P, freq)
        }
        P_matrix = rbind(P_matrix, P)
    }
    pbar1 = matrix(NA, nrow = n_cages, ncol = n_sites)
    pbar2 = matrix(NA, nrow = n_cages, ncol = n_sites)
    for (i in 1:n_cages) {
        pbar1[i, ] = P_matrix[((i - 1) * n_sites + 1):((i - 1) * 
            n_sites + n_sites), (ngen1 + 1)]
    }
    for (i in 1:n_cages) {
        pbar2[i, ] = P_matrix[((i - 1) * n_sites + 1):((i - 1) * 
            n_sites + n_sites), (ngen2 + 1)]
    }
    if (randomise) {
        if (verbose) {
            message("Randomising reference alleles in c_ind...")
        }
        ran_vect = sample(c(0, -1), ncol(c_genome), replace = T)
        ran_matrix = t(matrix(ran_vect, nrow = ncol(c_genome), 
            ncol = nrow(c_genome)))
        c_genome = abs(c_genome + ran_matrix)
        if (verbose) {
            message("Randomising reference alleles in list_alpha...")
        }
        ran_vect_alpha = ifelse(ran_vect == -1, -1, 1)
        list_alpha = list_alpha * ran_vect_alpha
        if (verbose) {
            message("Randomising reference alleles in pbar1 and pbar2...")
        }
        ran_matrix_pbar = t(matrix(ran_vect, nrow = ncol(pbar1), 
            ncol = nrow(pbar1)))
        pbar0 = abs(pbar0 + ran_vect)
        pbar1 = abs(pbar1 + ran_matrix_pbar)
        pbar2 = abs(pbar2 + ran_matrix_pbar)
    }
    return(list(c_genome = c_genome, list_alpha = list_alpha, 
        SNPs = SNPs, ngen1 = ngen1, ngen2 = ngen2, pbar0 = pbar0, 
        pbar1 = pbar1, pbar2 = pbar2, sim_params = sim_params[, 
            1:20]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
