\name{plot_alpha_distribution}
\alias{plot_alpha_distribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot properties of the joint distribution of alphas and allele frequencies
}
\description{
Plot properties of the joint distribution of alphas and allele frequencies
}
\usage{
plot_alpha_distribution(alpha = NULL, p = NULL, tprojp = NULL, parameters, ...)
}
\arguments{
  \item{alpha}{vector of average effects}
  \item{p}{vector of allele frequencies}
  \item{tprojp}{optional projection matrix (transposed) for allele frequencies and average effects}
  \item{parameters}{a list containing the elements \code{logit}, \code{balpha_0}, \code{balpha_1} and \code{palpha}. \code{logit} should be logical indicating whether the model is on the logit scale or not, and the remaining elements should be the numeric values of the parameters: see \link{}}
  \item{\dots}{additional arguments passed to the plot function}
}
\details{
}
\value{
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
function (alpha = NULL, p = NULL, tprojp = NULL, parameters, 
    ...) 
{
    pmq <- 2 * p - 1
    pq <- p * (1 - p)/2
    if (parameters$logit) {
        pred <- pmq_trans(pmq, parameters$balpha_0)
        intercept <- 0
    }
    else {
        pred <- pmq
        intercept <- parameters$balpha_0
    }
    slope <- parameters$balpha_1
    palpha <- parameters$palpha
    var <- pq^palpha
    if (is.null(tprojp)) {
        pmq.x <- seq(min(pmq), max(pmq), length = 1000)
        pq.x = (pmq.x + 1) * (1 - pmq.x)/8
        sd.x <- parameters$sigma2alpha * sqrt(pq.x^palpha)
        plot(alpha ~ pmq, ylab = "alpha", xlab = "p-q", ...)
        if (parameters$logit) {
            pred.x <- pmq_trans(pmq.x, parameters$balpha_0)
        }
        else {
            pred.x <- pmq.x
        }
        m.x <- intercept + slope * pred.x
        u.x <- qnorm(0.975, intercept + slope * pred.x, sd.x)
        l.x <- qnorm(0.025, intercept + slope * pred.x, sd.x)
        lines(m.x ~ pmq.x, col = "red")
        lines(u.x ~ pmq.x, col = "red", lty = 2)
        lines(l.x ~ pmq.x, col = "red", lty = 2)
    }
    else {
        alpha <- t(alpha \%*\% tprojp)
        pred <- t(pred \%*\% tprojp)
        var <- diag(t(tprojp) \%*\% (tprojp * pq^palpha))
        plot(alpha ~ pred, ylab = "projected alpha", xlab = "projected predictor", 
            ...)
        pred.x <- seq(min(pred), max(pred), length = 1000)
        m.x <- intercept + slope * pred.x
        lines(m.x ~ pred.x, col = "red")
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
