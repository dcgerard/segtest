% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polymapr.R
\name{polymapr_test}
\alias{polymapr_test}
\title{Run segregation distortion tests as implemented in the polymapR package.}
\usage{
polymapr_test(x, g1 = NULL, g2 = NULL, type = c("segtest", "polymapR"))
}
\arguments{
\item{x}{Either a vector of genotype counts, or a matrix of genotype
posteriors where the rows index the individuals and the columns
index the genotypes.}

\item{g1}{Parent 1's genotype.}

\item{g2}{Parent 2's genotype.}

\item{type}{Either my implementation which approximates that of
polymapR (\code{"segtest"}) or the implementation
through polymapR (\code{"polymapR"}). Note that
polymapR needs to be installed for \code{type = "polymapR"}.}
}
\value{
A list with the following elements:
\describe{
\item{p_value}{The p-value of the test.}
\item{bestfit}{The genotype frequencies of the best fit model.}
\item{frq_invalid}{The frequency of invalid genotypes.}
\item{p_invalid}{The p-value of the invalid proportion.}
}
}
\description{
The polymapR package tests for segregation distortion by iterating through all
possible forms of disomic or polysomic inheritance from either parent,
tests for concordance of the offspring genotypes using a chi-squared
test, and returns the largest p-value. It sometimes chooses a different
p-value based on other heuristics. They also sometimes return NA.
When \code{type = "segtest"}, we only look at patterns of the
given parent genotypes, choosing the largest p-value. When
\code{type = "polymapR"}, we return what they use via their heuristics.
}
\examples{
g1 <- 0
g2 <- 1
x <- c(4, 16, 0, 0, 0)
polymapr_test(x = x, g1 = g1, g2 = g2, type = "segtest")
polymapr_test(x = x, g1 = g1, g2 = g2, type = "polymapR")


}
\seealso{
\code{\link[polymapR]{checkF1}()}.
}
\author{
David Gerard
}
