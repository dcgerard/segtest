% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chisq.R
\name{chisq_g}
\alias{chisq_g}
\alias{chisq_g4}
\title{Chi Square test when genotypes are known}
\usage{
chisq_g(x, g1, g2)

chisq_g4(x, g1, g2)
}
\arguments{
\item{x}{Vector of observed genotype counts}

\item{g1}{Parent 1's genotype}

\item{g2}{Parent 2's genotype}
}
\value{
A list containing the chi-squared statistic, degrees of
freedom, and p-value.
}
\description{
This chi-squared test is run under the assumption of no double reduction
and no preferential pairing.
}
\section{Functions}{
\itemize{
\item \code{chisq_g4()}: Alias for chisq_g, for backwards compatibility.

}}
\examples{
x <- c(1, 2, 4, 3, 0)
g1 <- 2
g2 <- 2
chisq_g(x, g1, g2)

x <- c(10, 25, 10, 0, 0)
g1 <- 1
g2 <- 1
chisq_g(x, g1, g2)

}
\author{
Mira Thakkar and David Gerard
}
