
#' Parallelized likelihood ratio test for segregation distortion.
#'
#' @param g One of two inputs
#'   \itemize{
#'     \item{A matrix of genotype counts. The rows index the locis and the columns index the genotypes.}
#'     \item{An array of genotype log-likelihoods. The rows index the loci, the columns index the individuals, and the slices index the genotypes. Log-likelihoods are base e (natural log).}
#'   }
#' @param p1 One of two inputs
#'   \itemize{
#'     \item{A vector of parent 1's genotypes.}
#'     \item{A matrix of parent 1's genotype log-likelihoods. The rows index the loci and the columns index the genotypes. Logs are in base e (natural log).}
#'   }
#' @param p2 One of two inputs
#'   \itemize{
#'     \item{A vector of parent 1's genotypes.}
#'     \item{A matrix of parent 1's genotype log-likelihoods. The rows index the loci and the columns index the genotypes. Logs are in base e (natural log).}
#'   }
#' @inheritParams lrt_men_gl4
#'
#' @author David Gerard
#'
#' @export
multi_lrt <- function(g,
                      p1,
                      p2,
                      drbound = 1/6,
                      pp = TRUE,
                      dr = TRUE,
                      alpha = 0,
                      xi1 = 1/3,
                      xi2 = 1/3) {

}
