% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lrt_glpknown.R
\name{lrt_men_gl4}
\alias{lrt_men_gl4}
\title{Likelihood ratio test using genotype likelihoods.}
\usage{
lrt_men_gl4(
  gl,
  g1 = NULL,
  g2 = NULL,
  drbound = 1/6,
  pp = TRUE,
  dr = TRUE,
  alpha = 0,
  xi1 = 1/3,
  xi2 = 1/3
)
}
\arguments{
\item{gl}{The genotype log-likelihoods. The rows index the individuals
and the columns index the genotypes.}

\item{g1}{Either parent 1's genotype, or parent 1's genotype log-likelihoods.}

\item{g2}{Either parent 2's genotype, or parent 2's genotype log-likelihoods.}

\item{drbound}{The upper bound on the double reduction rate.}

\item{pp}{Is (partial) preferential pairing possible (\code{TRUE}) or
not (\code{FALSE})?}

\item{dr}{Is double reduction possible (\code{TRUE}) or
not (\code{FALSE})?}

\item{alpha}{If \code{dr = FALSE}, this is the known rate of double
reduction.}

\item{xi1}{If \code{pp = FALSE}, this is the known preferential pairing
parameter of parent 1.}

\item{xi2}{If \code{pp = FALSE}, this is the known preferential pairing
parameter of parent 2.}
}
\value{
A list with the following elements
\describe{
\item{\code{statistic}}{The log-likelihood ratio test statistic.}
\item{\code{df}}{The degrees of freedom.}
\item{\code{p_value}}{The p-value.}
\item{\code{alpha}}{The estimated double reduction rate.}
\item{\code{xi1}}{The estimated preferential pairing parameter of parent 1.}
\item{\code{xi2}}{The estimated preferential pairing parameter of parent 2.}
}
}
\description{
This will run a likelihood ratio test using the genotypes of an F1 population
of tetraploids for the null of Mendelian segregation (accounting for double
reduction and preferential pairing) against the alternative of
segregation distortion. This is when genotype uncertainty is accounted
for through genotype likelihoods.
}
\section{Unidentified parameters}{

When \code{g1 = 2} or \code{g2 = 2} (or both), the model is not identified
and those estimates (\code{alpha}, \code{xi1}, and \code{xi2}) are
meaningless. Do NOT interpret them.

The estimate of \code{alpha} (double reduction rate) IS identified as
long as at least one parent is simplex, and no parent is duplex.
However, the estimates of the double reduction rate have extremely high
variance.
}

\examples{
## null simulation
set.seed(1)
g1 <- 2
g2 <- 2
gl <- simf1gl(n = 25, g1 = g1, g2 = g2, alpha = 1/12, xi2 = 1/2)

## LRT when parent genotypes are known.
lrt_men_gl4(gl = gl, g1 = g1, g2 = g2)

## LRT when parent genotypes are not known
lrt_men_gl4(gl = gl)

## Alternative simulation
gl <- simgl(nvec = rep(5, 5))
lrt_men_gl4(gl = gl, g1 = g1, g2 = g2)

}
\author{
David Gerard
}
