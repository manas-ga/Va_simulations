\name{fit.model}
\alias{fit.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fit.model(palpha, balpha, LDalpha, nsnps, UL, DL, L, ngen2, ngen1, nind, tprojp, pbar0, pbar1, pbar2, nrep, Selec, Ne_factor = 1, LLonly = FALSE, method = "REML", verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{palpha}{
%%     ~~Describe \code{palpha} here~~
}
  \item{balpha}{
%%     ~~Describe \code{balpha} here~~
}
  \item{LDalpha}{
%%     ~~Describe \code{LDalpha} here~~
}
  \item{nsnps}{
%%     ~~Describe \code{nsnps} here~~
}
  \item{UL}{
%%     ~~Describe \code{UL} here~~
}
  \item{DL}{
%%     ~~Describe \code{DL} here~~
}
  \item{L}{
%%     ~~Describe \code{L} here~~
}
  \item{ngen2}{
%%     ~~Describe \code{ngen2} here~~
}
  \item{ngen1}{
%%     ~~Describe \code{ngen1} here~~
}
  \item{nind}{
%%     ~~Describe \code{nind} here~~
}
  \item{tprojp}{
%%     ~~Describe \code{tprojp} here~~
}
  \item{pbar0}{
%%     ~~Describe \code{pbar0} here~~
}
  \item{pbar1}{
%%     ~~Describe \code{pbar1} here~~
}
  \item{pbar2}{
%%     ~~Describe \code{pbar2} here~~
}
  \item{nrep}{
%%     ~~Describe \code{nrep} here~~
}
  \item{Selec}{
%%     ~~Describe \code{Selec} here~~
}
  \item{Ne_factor}{
%%     ~~Describe \code{Ne_factor} here~~
}
  \item{LLonly}{
%%     ~~Describe \code{LLonly} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
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
function (palpha, balpha, LDalpha, nsnps, UL, DL, L, ngen2, ngen1, 
    nind, tprojp, pbar0, pbar1, pbar2, nrep, Selec, Ne_factor = 1, 
    LLonly = FALSE, method = "REML", verbose = TRUE) 
{
    if (verbose) {
        message("Computing the covariance structure of locus effects...")
        message("Computing covp over 1 generation...")
    }
    if (LDalpha) {
        if (palpha == 0) {
            covp <- diag(nsnps)
        }
        else {
            covp <- UL \%*\% diag(DL^(2 * palpha)) \%*\% t(UL)
        }
    }
    else {
        covp <- diag(L)^palpha
    }
    rm(list = c("UL", "DL"))
    gc(verbose = FALSE)
    pmq <- 2 * pbar0 - 1
    pmq <- Selec \%*\% pmq
    int <- Selec \%*\% rep(1, nsnps)
    if (verbose) {
        message("Computing SC...")
    }
    projpSelec <- t(tprojp) \%*\% Selec
    rm("Selec")
    gc(verbose = FALSE)
    if (LDalpha) {
        SC <- projpSelec \%*\% covp
        rm("covp")
        gc(verbose = FALSE)
        SC <- SC \%*\% t(projpSelec)
    }
    else {
        SC <- projpSelec \%*\% (t(projpSelec) * covp)
    }
    attr(SC, "INVERSE") <- FALSE
    dimnames(SC) <- list(1:nrow(SC), 1:nrow(SC))
    if (verbose) {
        message("Projecting allele frequencies...")
    }
    pbar1_proj <- pbar1 \%*\% tprojp
    pbar2_proj <- pbar2 \%*\% tprojp
    pmq_proj <- t(matrix(pmq, nsnps, nrep)) \%*\% tprojp
    int_proj <- t(matrix(int, nsnps, nrep)) \%*\% tprojp
    dat.gaussian <- data.frame(delta = c(pbar2_proj - pbar1_proj), 
        locus = gl(ncol(pbar1_proj), nrep, ncol(pbar1_proj) * 
            nrep), rep = gl(nrep, 1, ncol(pbar1_proj) * nrep), 
        pmq = c(pmq_proj), int = c(int_proj))
    if (!is.na(balpha[1])) {
        dat.gaussian$int <- dat.gaussian$int * balpha[1]
    }
    if (!is.na(balpha[2])) {
        dat.gaussian$pmq <- dat.gaussian$pmq * balpha[2]
    }
    if (verbose) {
        message("Fitting linear mixed model...")
    }
    prior <- list(R = list(V = 1, nu = 0), G = list(G1 = list(V = 1, 
        nu = 1, alpha.mu = 0, alpha.V = 1000)))
    SC[lower.tri(SC)] <- t(SC)[lower.tri(SC)]
    if (method == "REML") {
        if (is.na(balpha[1]) & is.na(balpha[2])) {
            m1 <- asreml(delta ~ int + pmq - 1, random = ~vm(locus, 
                SC, singG = "PSD"), data = dat.gaussian, Cfixed = TRUE)
        }
        if (is.na(balpha[1]) & !is.na(balpha[2])) {
            m1 <- asreml(delta ~ int + offset(pmq) - 1, random = ~vm(locus, 
                SC, singG = "PSD"), data = dat.gaussian, Cfixed = TRUE)
        }
        if (!is.na(balpha[1]) & is.na(balpha[2])) {
            m1 <- asreml(delta ~ offset(int) + pmq - 1, random = ~vm(locus, 
                SC, singG = "PSD"), data = dat.gaussian, Cfixed = TRUE)
        }
        if (!is.na(balpha[1]) & !is.na(balpha[2])) {
            warning("asreml doesn't allow models without fixed effects, so intercept fitted but replaced with balpha[1]!")
            m1 <- asreml(delta ~ offset(pmq + int), random = ~vm(locus, 
                SC, singG = "PSD"), data = dat.gaussian, Cfixed = TRUE)
        }
    }
    if (method == "MCMC") {
        if (LLonly) {
            stop("method = MCMC specified so can't return log-likelihood with LLony=TRUE")
        }
        invSC <- solve(SC)
        invSC <- as(invSC, "sparseMatrix")
        attr(invSC, "rowNames") <- 1:retain
        attr(invSC, "colNames") <- 1:retain
        prior <- list(B = list(mu = balpha, V = diag(2) * 1e+10))
        diag(prior$B$V)[which(!is.na(balpha))] <- 1e-10
        prior$B$mu[which(is.na(balpha))] <- 0
        m1 <- MCMCglmm(delta ~ pmq + int - 1, random = ~locus, 
            data = dat.gaussian, ginverse = list(locus = invSC), 
            family = "gaussian", pr = TRUE, prior = prior)
    }
    if (LLonly) {
        return(m1$loglik)
    }
    else {
        return(list(data = dat.gaussian, model = m1, SC = SC, 
            palpha = palpha))
    }
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
