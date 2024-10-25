\name{analyse_parents}
\alias{analyse_parents}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
\ Obtain useful information about the parents' generation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
analyse_parents(c_genome, list_alpha, compute_svdL = FALSE, LDalpha = FALSE, SNPs, RecombRate, HapLength, tol = sqrt(.Machine$double.eps), AtleastOneRecomb = FALSE, calc_nR = TRUE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{c_genome}{
%%     ~~Describe \code{c_genome} here~~
}
  \item{list_alpha}{
%%     ~~Describe \code{list_alpha} here~~
}
  \item{compute_svdL}{
%%     ~~Describe \code{compute_svdL} here~~
}
  \item{LDalpha}{
%%     ~~Describe \code{LDalpha} here~~
}
  \item{SNPs}{
%%     ~~Describe \code{SNPs} here~~
}
  \item{RecombRate}{
%%     ~~Describe \code{RecombRate} here~~
}
  \item{HapLength}{
%%     ~~Describe \code{HapLength} here~~
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
}
  \item{AtleastOneRecomb}{
%%     ~~Describe \code{AtleastOneRecomb} here~~
}
  \item{calc_nR}{
%%     ~~Describe \code{calc_nR} here~~
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
function (c_genome, list_alpha, compute_svdL = FALSE, LDalpha = FALSE, 
    SNPs, RecombRate, HapLength, tol = sqrt(.Machine$double.eps), 
    AtleastOneRecomb = FALSE, calc_nR = TRUE, verbose = TRUE) 
{
    n0_individuals = nrow(c_genome)/2
    if (verbose) {
        message("Calculating L...")
    }
    paternal <- seq(1, 2 * n0_individuals, 2)
    maternal <- paternal + 1
    c0 <- (c_genome[paternal, ] + c_genome[maternal, ])/2
    L <- cov(c0) * (n0_individuals - 1)/n0_individuals
    if (compute_svdL | LDalpha) {
        if (verbose) {
            message("Performing SVD on c0...")
        }
        svdC <- svd(scale(sqrt(1/n0_individuals) * c0, scale = FALSE), 
            nu = 0)
        retain <- sum(svdC$d > tol)
        UL <- svdC$v[, 1:retain]
        DL <- svdC$d[1:retain]
    }
    if (calc_nR) {
        if (verbose) {
            message("Calculating the matrix of non-recombinant fractions...")
        }
        nR = form_nR(SNPs, RecombRate, HapLength, AtleastOneRecomb)
        Lgp <- (cov(c_genome[paternal, ]) + cov(c_genome[maternal, 
            ])) * (n0_individuals - 1)/(4 * n0_individuals)
        Ltilde <- Lgp + (1 - nR) * (L - Lgp)/nR
        rm("Lgp")
    }
    else {
        nR = NULL
        Ltilde = NULL
    }
    seg_sites = length(list_alpha)
    seg_sites_neu = sum(list_alpha == 0)
    seg_sites_ben = sum(list_alpha > 0)
    seg_sites_del = sum(list_alpha < 0)
    pbar0 = colMeans(c_genome)
    diversity = pbar0 * (1 - pbar0)/2
    mean_diversity = mean(diversity)
    if (verbose) {
        message("Calculating the true additive genetic and genic variances...")
    }
    va_true = sum(diversity * list_alpha^2)
    vA_true = t(list_alpha) \%*\% L \%*\% list_alpha
    if (verbose) {
        message("Calculating properties of the empirical distribution of alphas...")
    }
    alpha_properties = alpha_distribution(alpha = list_alpha, 
        p = pbar0)
    palpha_emp = alpha_properties$palpha
    balpha_intercept_emp = alpha_properties$balpha_0
    balpha_slope_emp = alpha_properties$balpha_1
    sigma2alpha_emp = alpha_properties$sigma2alpha
    balpha_emp = c(balpha_intercept_emp, balpha_slope_emp)
    if (LDalpha) {
        TrV <- sum(DL^(2 * (palpha_emp + 1))) * sigma2alpha_emp
        aLa <- t((alpha_properties$X) \%*\% balpha_emp) \%*\% L \%*\% 
            (alpha_properties$X) \%*\% balpha_emp - sum(diag(t(alpha_properties$X) \%*\% 
            L \%*\% (alpha_properties$X) \%*\% (alpha_properties$S)))
    }
    else {
        TrV <- sum(diag(L)^(palpha_emp + 1)) * sigma2alpha_emp
        aLa <- t((alpha_properties$X) \%*\% balpha_emp) \%*\% L \%*\% 
            (alpha_properties$X) \%*\% balpha_emp - sum(diag(t(alpha_properties$X) \%*\% 
            L \%*\% (alpha_properties$X) \%*\% (alpha_properties$S)))
    }
    vA_alpha_emp <- TrV + aLa
    return(list(L = L, Ltilde = Ltilde, nR = nR, svdL = if (compute_svdL | 
        LDalpha) {
        list(UL = UL, DL = DL)
    } else {
        NULL
    }, seg_sites = seg_sites, seg_sites_neu = seg_sites_neu, 
        seg_sites_ben = seg_sites_ben, seg_sites_del = seg_sites_del, 
        mean_diversity = mean_diversity, va_true = va_true, vA_true = vA_true, 
        vA_alpha_emp = vA_alpha_emp, parameters = alpha_properties[1:5]))
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
