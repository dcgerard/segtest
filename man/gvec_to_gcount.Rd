% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rscripts.R
\name{gvec_to_gcount}
\alias{gvec_to_gcount}
\title{Inverse function of \code{\link{gcount_to_gvec}()}.}
\usage{
gvec_to_gcount(gvec, ploidy = 4)
}
\arguments{
\item{gvec}{The vector of genotypes. \code{gvec[i]} is the genotype
for individual i.}

\item{ploidy}{The ploidy of the species.}
}
\description{
Inverse function of \code{\link{gcount_to_gvec}()}.
}
\examples{
gvec <- c(1, 2, 3, 2, 3, 1, 4, 0, 1, 0, 0, 1, 0, 0)
gvec_to_gcount(gvec = gvec)

}
\seealso{
\code{\link[=gcount_to_gvec]{gcount_to_gvec()}}
}
\author{
David Gerard
}
