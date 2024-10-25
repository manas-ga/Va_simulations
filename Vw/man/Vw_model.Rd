\name{Vw_model}
\alias{Vw_model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Vw_model(c_genome, nR, pbar0, pbar1, ngen1 = 1, pbar2, ngen2, nind, proj, LDalpha, pa, palpha, balpha, Vs, method, L, Ltilde, svdL = NULL, Ne_factor = 1, tol = sqrt(.Machine$double.eps), save_tprojp = FALSE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{c_genome}{
%%     ~~Describe \code{c_genome} here~~
}
  \item{nR}{
%%     ~~Describe \code{nR} here~~
}
  \item{pbar0}{
%%     ~~Describe \code{pbar0} here~~
}
  \item{pbar1}{
%%     ~~Describe \code{pbar1} here~~
}
  \item{ngen1}{
%%     ~~Describe \code{ngen1} here~~
}
  \item{pbar2}{
%%     ~~Describe \code{pbar2} here~~
}
  \item{ngen2}{
%%     ~~Describe \code{ngen2} here~~
}
  \item{nind}{
%%     ~~Describe \code{nind} here~~
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
  \item{palpha}{
%%     ~~Describe \code{palpha} here~~
}
  \item{balpha}{
%%     ~~Describe \code{balpha} here~~
}
  \item{Vs}{
%%     ~~Describe \code{Vs} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{L}{
%%     ~~Describe \code{L} here~~
}
  \item{Ltilde}{
%%     ~~Describe \code{Ltilde} here~~
}
  \item{svdL}{
%%     ~~Describe \code{svdL} here~~
}
  \item{Ne_factor}{
%%     ~~Describe \code{Ne_factor} here~~
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
}
  \item{save_tprojp}{
%%     ~~Describe \code{save_tprojp} here~~
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
function (c_genome, nR, pbar0, pbar1, ngen1 = 1, pbar2, ngen2, 
    nind, proj, LDalpha, pa, palpha, balpha, Vs, method, L, Ltilde, 
    svdL = NULL, Ne_factor = 1, tol = sqrt(.Machine$double.eps), 
    save_tprojp = FALSE, verbose = TRUE) 
{
    asreml.options(Cfixed = TRUE)
    if (is.null(L) | is.null(Ltilde) | (is.null(svdL) & (proj == 
        "BLoM" | LDalpha))) {
        if (is.null(c_genome)) {
            stop("c_genome is required if L, Ltilde or svdL are NULL")
        }
        n0_individuals <- nrow(c_genome)/2
        if (is.null(L) | (is.null(svdL) & (proj == "BLoM" | LDalpha))) {
            paternal <- seq(1, 2 * n0_individuals, 2)
            maternal <- paternal + 1
            c0 <- (c_genome[paternal, ] + c_genome[maternal, 
                ])/2
            if (is.null(L)) {
                if (verbose) {
                  message("Computing L in the parents' generation...")
                }
                L <- cov(c0) * (n0_individuals - 1)/n0_individuals
            }
        }
        if (is.null(Ltilde)) {
            Lgp <- (cov(c_genome[paternal, ]) + cov(c_genome[maternal, 
                ])) * (n0_individuals - 1)/(4 * n0_individuals)
            Ltilde <- Lgp + (1 - nR) * (L - Lgp)/nR
            rm("Lgp")
        }
    }
    if (!proj \%in\% c("LoM", "BLoM", "L", "N")) {
        stop("proj must be one of 'LoM', 'L', 'N'")
    }
    if (!Vs \%in\% c("LoNL", "L")) {
        stop("Vs must be either 'LoNL' or 'L'")
    }
    if (!method \%in\% c("REML", "MCMC")) {
        stop("method must be either 'REML' or 'MCMC'")
    }
    nsnps <- ncol(pbar1)
    nrep <- nrow(pbar1)
    if (ncol(pbar1) != nsnps) {
        stop("pbar1 should have as many columns as the number of columns in Co (the number of SNPs)")
    }
    if (ncol(pbar2) != nsnps) {
        stop("pbar2 should have as many columns as the number of columns in Co (the number of SNPs)")
    }
    if (nrow(pbar1) != nrow(pbar2)) {
        stop("pbar1 and pbar2 should have the same number of rows (the number of replicates)")
    }
    if (proj == "BLoM" | LDalpha) {
        if (is.null(svdL)) {
            if (verbose) {
                message("Performing SVD on C0...")
            }
            c0 = scale(sqrt(1/nind) * c0, scale = FALSE)
            svdC <- svd(c0, nu = 0)
            retain <- sum(svdC$d > tol)
            UL <- svdC$v[, 1:retain]
            DL <- svdC$d[1:retain]
            rm("svdC")
            gc(verbose = FALSE)
        }
        else {
            UL <- svdL$UL
            DL <- svdL$DL
            rm("svdL")
            gc(verbose = FALSE)
            retain <- ncol(UL)
        }
    }
    if (proj == "LoM" | proj == "BLoM") {
        if (verbose) {
            message("Computing M...")
        }
        M = ((1 - 1/(2 * nind * Ne_factor))^(ngen1)) * (1/(nind * 
            Ne_factor)) * nR^(1 + ngen1)
        if ((ngen2 - ngen1) > 1) {
            for (x in (ngen1 + 2):(ngen2)) {
                M = M + ((1 - 1/(2 * nind * Ne_factor))^(x - 
                  1)) * (1/(nind * Ne_factor)) * nR^x
            }
        }
        Drift <- Ltilde * M
        rm("M")
        gc(verbose = FALSE)
    }
    if (Vs == "LoNL") {
        if (verbose) {
            message("Computing N...")
        }
        if (ngen1 != 0) {
            N = {
                ((1 - 1/(2 * nind * Ne_factor))^(ngen1)) * nR^(ngen1)
            }
        }
        else {
            if (ngen2 > 1) {
                N = matrix(0, ncol(nR), ncol(nR))
            }
        }
        if ((ngen2 - ngen1) > 1) {
            for (x in (ngen1 + 2):(ngen2)) {
                N = N + {
                  ((1 - 1/(2 * nind * Ne_factor))^(x - 1)) * 
                    nR^(x - 1)
                }
            }
        }
        if (ngen1 == 0) {
            Selec <- L
            if (ngen2 > 1) {
                Selec <- Selec + Ltilde * N
                rm("N")
                gc(verbose = FALSE)
            }
        }
        else {
            Selec <- Ltilde * N
            rm("N")
            gc(verbose = FALSE)
        }
    }
    else {
        Selec <- L * (ngen2 - ngen1)
    }
    if (proj == "LoM") {
        if (verbose) {
            message("Performing eigendecomposition of L*M...")
        }
        sdLoM <- RSpectra::eigs(Drift, min(nind, nsnps))
        retain <- sum(sqrt(sdLoM$values) > tol)
        U <- sdLoM$vectors[, 1:retain]
        D <- sqrt(sdLoM$values[1:retain])
        rm("sdLoM")
        gc(verbose = FALSE)
    }
    if (proj == "BLoM") {
        if (verbose) {
            message("Performing eigendecomposition of t(U)\%*\%(L*M)\%*\%U...")
        }
        sdLoM <- RSpectra::eigs(t(UL) \%*\% Drift \%*\% UL, retain)
        retain <- sum(sqrt(sdLoM$values) > tol)
        U <- UL \%*\% sdLoM$vectors[, 1:retain]
        D <- sqrt(sdLoM$values[1:retain])
        rm("sdLoM")
        gc(verbose = FALSE)
    }
    if (verbose) {
        message("Computing the projection matrix...")
    }
    if (proj == "N") {
        tprojp <- diag(nrow(pbar1))
    }
    else {
        tprojp <- U \%*\% diag(D^(-pa))
    }
    rm(list = c("U", "D"))
    gc(verbose = FALSE)
    if (is.na(palpha)) {
        if (verbose) {
            message("Estimating palpha...")
        }
        palpha <- optim(0, fit.model, balpha = balpha, LDalpha = LDalpha, 
            nsnps = nsnps, UL = UL, DL = DL, L = L, ngen2 = ngen2, 
            ngen1 = ngen1, nind = nind, tprojp = tprojp, pbar0 = pbar0, 
            pbar1 = pbar1, pbar2 = pbar2, nrep = nrep, LLonly = TRUE, 
            Selec = Selec, Ne_factor = Ne_factor, verbose = verbose, 
            method = "L-BFGS-B", lower = -2, upper = 2, control = list(fnscale = -1, 
                factr = 1e+11), hessian = TRUE)
        palpha_var <- -1/palpha$hessian
        palpha <- palpha$par
    }
    else {
        palpha_var <- 0
    }
    if (verbose) {
        message("Fitting the final model...")
    }
    output <- fit.model(palpha = palpha, balpha = balpha, LDalpha = LDalpha, 
        nsnps = nsnps, UL = UL, DL = DL, L = L, ngen2 = ngen2, 
        ngen1 = ngen1, nind = nind, tprojp = tprojp, pbar0 = pbar0, 
        pbar1 = pbar1, pbar2 = pbar2, nrep = nrep, LLonly = FALSE, 
        Selec = Selec, Ne_factor = Ne_factor, verbose = verbose)
    if (verbose) {
        message("Calculating the estimate of Vw...")
    }
    if (method == "REML") {
        sigma2alpha <- summary(output$model)$varcomp[1, 1]
        S <- matrix(0, 2, 2)
        if (is.na(balpha[1]) & is.na(balpha[2])) {
            balpha <- summary(output$model, coef = TRUE)$coef.fixed[c("int", 
                "pmq"), 1]
            S <- output$model$Cfixed[c("int", "pmq"), c("int", 
                "pmq")]
        }
        if (is.na(balpha[1]) & !is.na(balpha[2])) {
            balpha[1] <- summary(output$model, coef = TRUE)$coef.fixed["int", 
                1]
            S["int", "int"] <- output$model$Cfixed
        }
        if (!is.na(balpha[1]) & is.na(balpha[2])) {
            balpha[2] <- summary(output$model, coef = TRUE)$coef.fixed["pmq", 
                1]
            S["pmq", "pmq"] <- output$model$Cfixed
        }
    }
    if (method == "MCMC") {
        sigma2alpha <- posterior.mode(output$model$VCV[, 1])
        balpha <- colMeans(output$model$Sol[, c("int", "pmq")])
    }
    X <- cbind(rep(1, length(pbar0)), 2 * pbar0 - 1)
    if (LDalpha) {
        TrV <- sum(DL^(2 * (output$palpha + 1))) * sigma2alpha
        aLa <- t(X \%*\% balpha) \%*\% L \%*\% X \%*\% balpha - sum(diag(t(X) \%*\% 
            L \%*\% X \%*\% S))
    }
    else {
        TrV <- sum(diag(L)^(output$palpha + 1)) * sigma2alpha
        aLa <- t(X \%*\% balpha) \%*\% L \%*\% X \%*\% balpha - sum(diag(t(X) \%*\% 
            L \%*\% X \%*\% S))
    }
    Vw_est <- TrV + aLa
    return(list(Vw_est = Vw_est, data = output$data, model = output$model, 
        SC = output$SC, palpha = output$palpha, balpha = balpha, 
        palpha_var = palpha_var, balpha_var = S, tprojp = if (save_tprojp) {
            tprojp
        } else {
            NULL
        }, S = S, X = X, DL = ifelse(exists("DL"), DL, NA)))
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
