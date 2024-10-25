\name{sim_exp}
\alias{sim_exp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
\sim exp
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sim_exp(c_genome, nind, r, position, nrep, alpha = NULL, Ve = 0, fitness_model = "Exp")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{c_genome}{
%%     ~~Describe \code{c_genome} here~~
}
  \item{nind}{
%%     ~~Describe \code{nind} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{position}{
%%     ~~Describe \code{position} here~~
}
  \item{nrep}{
%%     ~~Describe \code{nrep} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{Ve}{
%%     ~~Describe \code{Ve} here~~
}
  \item{fitness_model}{
%%     ~~Describe \code{fitness_model} here~~
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
function (c_genome, nind, r, position, nrep, alpha = NULL, Ve = 0, 
    fitness_model = "Exp") 
{
    if (ncol(c_genome) != length(position)) {
        stop("position should be the same length as number of columns in c_genome")
    }
    if (!is.null(alpha)) {
        if (ncol(c_genome) != length(alpha)) {
            stop("alpha should be the same length as number of columns in c_genome")
        }
    }
    reordering <- order(position)
    c_genome <- c_genome[, reordering]
    position <- position[reordering]
    n0_individuals <- nrow(c_genome)/2
    n_loci <- ncol(c_genome)
    pbar <- matrix(NA, nrep, n_loci)
    if (!is.null(alpha)) {
        alpha <- alpha[reordering]
        genetic_fitness <- c_genome \%*\% alpha/2
        genetic_fitness <- genetic_fitness[seq(1, 2 * n0_individuals, 
            2)] + genetic_fitness[seq(2, 2 * n0_individuals, 
            2)]
    }
    else {
        genetic_fitness <- rep(0, n0_individuals)
    }
    genetic_fitness <- 1 + genetic_fitness - mean(genetic_fitness)
    if (any(genetic_fitness < 0)) {
        warning(paste(sum(genetic_fitness < 0), "genetic fitnesses negative - set to zero"))
        genetic_fitness[which(genetic_fitness < 0)] <- 0
    }
    for (rep in 1:nrep) {
        fitness <- genetic_fitness
        if (Ve != 0) {
            if (fitness_model == "Exp") {
                fitness <- fitness + rnorm(n0_individuals, 0, 
                  sqrt(Ve))
            }
            if (fitness_model == "I") {
                m <- min(fitness)
                fitness <- fitness + rgamma(n0_individuals, shape = (m^2)/Ve, 
                  scale = Ve/m) - m
            }
        }
        if (fitness_model == "Exp") {
            fitness <- exp(fitness)
        }
        genomes <- matrix(0, nind, n_loci)
        for (i in 1:nind) {
            for (j in 1:2) {
                parent <- sample(1:n0_individuals, 1, prob = fitness)
                gamete <- rep(rbinom(1, prob = 0.5, size = 1), 
                  n_loci)
                r_event <- reda::simEventData(rho = r, end = max(position))
                if (nrow(r_event) > 0) {
                  r_event <- r_event$time[which(r_event$event == 
                    1)]
                  for (k in 1:length(r_event)) {
                    gamete[which(position > r_event[k])] <- abs(1 - 
                      gamete[which(position > r_event[k])])
                  }
                }
                gam0 <- which(gamete == 0)
                if (length(gam0) != 0) {
                  genomes[i, gam0] <- genomes[i, gam0] + c_genome[(parent - 
                    1) * 2 + 1, gam0]
                }
                if (length(gam0) != n_loci) {
                  if (length(gam0) != 0) {
                    genomes[i, -gam0] <- genomes[i, -gam0] + 
                      c_genome[(parent - 1) * 2 + 2, -gam0]
                  }
                  else {
                    genomes[i, ] <- genomes[i, ] + c_genome[(parent - 
                      1) * 2 + 2, ]
                  }
                }
            }
        }
        pbar[rep, ] <- colMeans(genomes)/2
    }
    pbar <- pbar[, order(reordering)]
    return(pbar)
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
