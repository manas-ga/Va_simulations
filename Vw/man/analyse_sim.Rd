\name{analyse_sim}
\alias{analyse_sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
\Analyse a SLiM simulation to estimate Vw using the method 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
analyse_sim(Set_ID, sim = 1, unzip = FALSE, slim_output_path, sim_param_path, extract_genomes_path, extract_mut_path, mutations_path, c_matrix_path, output_path, n_sample = NULL, randomise = TRUE, delete_temp_files = TRUE, proj = "BLoM", LDalpha = FALSE, pa = 1, Vs = "LoNL", method = "REML", palpha = NA, balpha = c(NA, NA), AtleastOneRecomb = FALSE, all.gp = FALSE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Set_ID}{
%%     ~~Describe \code{Set_ID} here~~
}
  \item{sim}{
%%     ~~Describe \code{sim} here~~
}
  \item{unzip}{
%%     ~~Describe \code{unzip} here~~
}
  \item{slim_output_path}{
%%     ~~Describe \code{slim_output_path} here~~
}
  \item{sim_param_path}{
%%     ~~Describe \code{sim_param_path} here~~
}
  \item{extract_genomes_path}{
%%     ~~Describe \code{extract_genomes_path} here~~
}
  \item{extract_mut_path}{
%%     ~~Describe \code{extract_mut_path} here~~
}
  \item{mutations_path}{
%%     ~~Describe \code{mutations_path} here~~
}
  \item{c_matrix_path}{
%%     ~~Describe \code{c_matrix_path} here~~
}
  \item{output_path}{
%%     ~~Describe \code{output_path} here~~
}
  \item{n_sample}{
%%     ~~Describe \code{n_sample} here~~
}
  \item{randomise}{
%%     ~~Describe \code{randomise} here~~
}
  \item{delete_temp_files}{
%%     ~~Describe \code{delete_temp_files} here~~
}
  \item{proj}{
%%     ~~Describe \code{proj} here~~
}
  \item{LDalpha}{
%%     ~~Describe \code{LDalpha} here~~
}
  \item{pa}{
%%     ~~Describe \code{pa} here~~
}
  \item{Vs}{
%%     ~~Describe \code{Vs} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{palpha}{
%%     ~~Describe \code{palpha} here~~
}
  \item{balpha}{
%%     ~~Describe \code{balpha} here~~
}
  \item{AtleastOneRecomb}{
%%     ~~Describe \code{AtleastOneRecomb} here~~
}
  \item{all.gp}{
%%     ~~Describe \code{all.gp} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
    output_path, n_sample = NULL, randomise = TRUE, delete_temp_files = TRUE, 
    proj = "BLoM", LDalpha = FALSE, pa = 1, Vs = "LoNL", method = "REML", 
    palpha = NA, balpha = c(NA, NA), AtleastOneRecomb = FALSE, 
    all.gp = FALSE, verbose = TRUE) 
{
    if (verbose) {
        message("Extracting simulation data...")
    }
    sim_data = extract_slim_data(Set_ID = Set_ID, sim = sim, 
        unzip = unzip, slim_output_path = slim_output_path, sim_param_path = sim_param_path, 
        extract_genomes_path = extract_genomes_path, extract_mut_path = extract_mut_path, 
        mutations_path = mutations_path, c_matrix_path = c_matrix_path, 
        randomise = randomise)
    if (verbose) {
        message("Analysing the parents' generation...")
    }
    parents_info = analyse_parents(c_genome = sim_data$c_genome, 
        list_alpha = sim_data$list_alpha, compute_svdL = TRUE, 
        LDalpha = LDalpha, SNPs = sim_data$SNPs, RecombRate = sim_data$sim_params$r_expt, 
        HapLength = sim_data$sim_params$sequence_length, AtleastOneRecomb = AtleastOneRecomb)
    if (verbose) {
        message("Performing analyses...")
    }
    m1 <- Vw_model(c_genome = NULL, nR = parents_info$nR, pbar0 = sim_data$pbar0, 
        pbar1 = sim_data$pbar1, ngen1 = sim_data$ngen1, pbar2 = sim_data$pbar2, 
        ngen2 = sim_data$ngen2, nind = sim_data$sim_params$n_ind_exp, 
        proj = proj, LDalpha = LDalpha, pa = pa, palpha = palpha, 
        balpha = balpha, Vs = Vs, method = method, L = parents_info$L, 
        Ltilde = if (all.gp) {
            parents_info$L
        }
        else {
            parents_info$Ltilde
        }, svdL = parents_info$svdL, tol = sqrt(.Machine$double.eps))
    vA_est = m1$Vw_est
    palpha_est = m1$palpha
    palpha_var_est = m1$palpha_var
    balpha_intercept_est = m1$balpha[1]
    balpha_slope_est = m1$balpha[2]
    balpha_var_est = paste(m1$balpha_var[1, 1], m1$balpha_var[2, 
        2], m1$balpha_var[1, 2], sep = "_")
    sigma2alpha_est = summary(m1$model)$varcomp[1, 1]
    if (verbose) {
        message("Saving data...")
    }
    sim_params = sim_data$sim_params
    analysis_data = data.frame(proj = proj, LDalpha = LDalpha, 
        pa = pa, Vs = Vs, randomise = randomise, palpha_method = palpha, 
        balpha_method = paste(balpha[1], balpha[2], sep = "_"), 
        va_true = parents_info$va_true, vA_true = parents_info$vA_true, 
        vA_est = vA_est, vA_alpha_emp = parents_info$vA_alpha_emp, 
        palpha_emp = parents_info$parameters$palpha, balpha_intercept_emp = parents_info$parameters$balpha_0, 
        balpha_slope_emp = parents_info$parameters$balpha_1, 
        sigma2alpha_emp = parents_info$parameters$sigma2alpha, 
        palpha_est = palpha_est, palpha_var_est = palpha_var_est, 
        balpha_intercept_est = balpha_intercept_est, balpha_slope_est = balpha_slope_est, 
        balpha_var_est = balpha_var_est, sigma2alpha_est = sigma2alpha_est, 
        seg_sites = parents_info$seg_sites, seg_sites_neu = parents_info$seg_sites_neu, 
        seg_sites_ben = parents_info$seg_sites_ben, seg_sites_del = parents_info$seg_sites_del, 
        mean_diversity = parents_info$mean_diversity, all.gp = all.gp)
    analysis_data = cbind(sim_params, analysis_data)
    write.table(rbind(names(analysis_data), analysis_data), file = paste(output_path, 
        "/", Set_ID, "_sim_", sim, "_Data_analysis.csv", sep = ""), 
        col.names = FALSE, row.names = FALSE, sep = ",")
    return(analysis_data)
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
