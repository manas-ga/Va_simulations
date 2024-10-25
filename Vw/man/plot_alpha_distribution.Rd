\name{plot_alpha_distribution}
\alias{plot_alpha_distribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_alpha_distribution(alpha = NULL, p = NULL, tprojp = NULL, parameters, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{tprojp}{
%%     ~~Describe \code{tprojp} here~~
}
  \item{parameters}{
%%     ~~Describe \code{parameters} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
