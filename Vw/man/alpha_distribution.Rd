\name{alpha_distribution}
\alias{alpha_distribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
\ Analyse properties of the joint distribution of alphas and genetic diversities
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
alpha_distribution(alpha, p, tprojp = NULL, logit = FALSE, save_model = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{
A vector of average effects for relative fitness
}
  \item{p}{
A vector of allele frequencies in the parents' generation
}
  \item{tprojp}{
%%     ~~Describe \code{tprojp} here~~
}
  \item{logit}{
%%     ~~Describe \code{logit} here~~
}
  \item{save_model}{
%%     ~~Describe \code{save_model} here~~
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
function (alpha, p, tprojp = NULL, logit = FALSE, save_model = FALSE) 
{
    pmq <- 2 * p - 1
    pq <- p * (1 - p)/2
    palpha_fit <- function(palpha, alpha, pmq, pq, tprojp) {
        if (!is.null(tprojp)) {
            alpha <- t(alpha \%*\% tprojp)
            pmq <- t(pmq \%*\% tprojp)
            pq <- diag(t(tprojp) \%*\% (tprojp * pq^palpha))
        }
        else {
            pq <- pq^palpha
        }
        -logLik(lm(alpha ~ pmq, weights = 1/pq))
    }
    logit_fit <- function(par, alpha, pmq, pq, tprojp) {
        balpha_0 <- par[1]
        palpha <- par[2]
        pmq <- pmq_trans(pmq, balpha_0)
        if (!is.null(tprojp)) {
            alpha <- t(alpha \%*\% tprojp)
            pmq <- t(pmq \%*\% tprojp)
            pq <- diag(t(tprojp) \%*\% (tprojp * pq^palpha))
        }
        else {
            pq <- pq^palpha
        }
        -logLik(lm(alpha ~ pmq - 1, weights = 1/pq))
    }
    if (logit) {
        palpha <- optim(c(1, 1), logit_fit, alpha = alpha, pmq = pmq, 
            pq = pq, tprojp = tprojp)$par
        balpha_0 <- palpha[1]
        palpha <- palpha[2]
        pmq <- pmq_trans(pmq, balpha_0)
        if (!is.null(tprojp)) {
            alpha <- t(alpha \%*\% tprojp)
            pmq <- t(pmq \%*\% tprojp)
            pq <- diag(t(tprojp) \%*\% (tprojp * pq^palpha))
        }
        else {
            pq <- pq^palpha
        }
        m1 <- lm(alpha ~ pmq - 1, weights = 1/pq)
    }
    else {
        palpha <- optim(0, palpha_fit, alpha = alpha, pmq = pmq, 
            pq = pq, tprojp = tprojp, method = "Brent", lower = -5, 
            upper = 5)$par
        if (!is.null(tprojp)) {
            alpha <- t(alpha \%*\% tprojp)
            pmq <- t(pmq \%*\% tprojp)
            pq <- diag(t(tprojp) \%*\% (tprojp * pq^palpha))
        }
        else {
            pq <- pq^palpha
        }
        m1 <- lm(alpha ~ pmq, weights = 1/pq)
    }
    result <- list(palpha = palpha, balpha_0 = if (logit) {
        balpha_0
    } else {
        coef(m1)["(Intercept)"]
    }, balpha_1 = coef(m1)["pmq"], sigma2alpha = summary(m1)$sigma^2, 
        logit = logit, model = if (save_model) {
            m1
        } else {
            NULL
        }, X = model.matrix(m1), S = summary(m1)$cov.unscaled * 
            summary(m1)$sigma^2)
    return(result)
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
