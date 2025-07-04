% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multi.R
\name{multidog_to_g}
\alias{multidog_to_g}
\title{Converts multidog output to a format usable for seg_multi() and multi_lrt()}
\usage{
multidog_to_g(
  mout,
  ploidy,
  type = c("off_gl", "all_gl", "off_g", "all_g"),
  p1 = NULL,
  p2 = NULL
)
}
\arguments{
\item{mout}{The output of \code{\link[updog]{multidog}()}.}

\item{ploidy}{The ploidy.}

\item{type}{\describe{
\item{\code{"off_gl"}}{Genotype likelihoods of offspring but not parents. This is the typical choice if you used the "f1", "f1pp", "s1", or "s1pp" options when genotyping.}
\item{\code{"all_gl"}}{Genotype likelihoods of offspring and parents. This is only done if you did \emph{not} use the "f1", "f1pp", "s1", or "s1pp" options when genotyping. If this is the case, then you need to specify which individuals are the parents.}
\item{\code{"off_g"}}{Genotypes, assuming that they are known. You used the "f1", "f1pp", "s1", or "s1pp" option when genotyping.}
\item{\code{"all_g"}}{Genotypes, assuming that they are known. You did \emph{not} use the "f1", "f1pp", "s1", or "s1pp" option when genotyping. If this is the case, then you need to specify which individuals are the parents.}
}}

\item{p1}{The first (or only) parent name if using \code{type = "all_gl"} or \code{type = "all_g"}.}

\item{p2}{The second parent name if using \code{type = "all_gl"} or \code{type = "all_g"}. Omit if you used the "s1" or "s1pp" models when genotyping.}
}
\value{
A list with the following elements
\describe{
\item{\code{g}}{Either a matrix of counts, where the columns index the genotype
and the rows index the loci (\code{type = "all_g"} or
\code{type = "off_g"}). Or an array of genotype (natural) log-likelihoods
where the rows index the loci, the columns index the
individuals, and the slices index the genotypes
(\code{type = "all_gl"} or \code{type = "off_gl"}).}
\item{\code{p1}}{Either a vector of known parental genotypes
(\code{type = "off_gl"}, \code{type = "all_g"} or \code{type = "off_g"}).
Or a matrix of genotype (natural) log-likelihoods where the
rows index the loci and the columns index the genotypes
(\code{type = "all_gl"}).}
\item{\code{p2}}{Either a vector of known parental genotypes
(\code{type = "off_gl"}, \code{type = "all_g"} or \code{type = "off_g"}).
Or a matrix of genotype (natural) log-likelihoods where the
rows index the loci and the columns index the genotypes
(\code{type = "all_gl"}).
This will be \code{NULL} if you (i) used \code{"s1"} or
\code{"s1pp"} models in updog and used either
\code{type = "off_g"} or \code{type = "off_gl"} or
(ii) used \code{type = "all_g"} or \code{type = "all_gl"}
and only specified \code{p1} but not \code{p2}.}
}
}
\description{
Converts multidog output to a format usable for seg_multi() and multi_lrt()
}
\examples{
multidog_to_g(
  mout = ufit,
  ploidy = 4,
  type = "all_g",
  p1 = "indigocrisp",
  p2 = "sweetcrisp")
multidog_to_g(
  mout = ufit,
  ploidy = 4,
  type = "all_gl",
  p1 = "indigocrisp",
  p2 = "sweetcrisp")
multidog_to_g(mout = ufit2, ploidy = 4, type = "off_g")
multidog_to_g(mout = ufit2, ploidy = 4, type = "off_gl")
multidog_to_g(mout = ufit3, ploidy = 4, type = "off_g")
multidog_to_g(mout = ufit3, ploidy = 4, type = "off_gl")


}
\author{
David Gerard
}
