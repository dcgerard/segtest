## Classic chi-squared tests

## Chi-Sq for GL
#' Chi-Sq for GL
#'
#' Calculates the MLE genotype and runs a chi-squared test assuming
#' no double reduction and no preferential pairing.
#'
#' @param gl A matrix of offspring genotype log-likelihoods. The rows index the
#'     individuals and the columns index the possible genotypes. So
#'     \code{gl[i, k]} is the offspring genotype log-likelihood for individual i and
#'     genotype k-1.
#' @param g1 The first parent's genotype.
#' @param g2 The second parent's genotype.
#'
#' @inherit chisq_g return
#'
#' @author Mira Thakkar and David Gerard
#'
#' @examples
#' ## null sim
#' set.seed(1)
#' g1 <- 2
#' g2 <- 2
#' gl <- simf1gl(n = 25, g1 = g1, g2 = g2, alpha = 0, xi2 = 1/3)
#' chisq_gl(gl = gl, g1 = g1, g2 = g2)
#'
#' @export
chisq_gl <- function(gl, g1, g2){
  ploidy <- ncol(gl) - 1
  col_max <- apply(gl, 1, which.max) - 1
  col_max <- factor(col_max, levels = 0:ploidy)
  x <- c(table(col_max))
  output <- chisq_g4(x = x, g1 = g1, g2 = g2)
  return(output)
}

#' Chi Square test when genotypes are known
#'
#' This chi-squared test is run under the assumption of no double reduction
#' and no preferential pairing.
#'
#' @param x Vector of observed genotype counts
#' @param g1 Parent 1's genotype
#' @param g2 Parent 2's genotype
#'
#' @return A list containing the chi-squared statistic, degrees of
#'     freedom, and p-value.
#'
#' @author Mira Thakkar and David Gerard
#'
#' @examples
#' x <- c(1, 2, 4, 3, 0)
#' g1 <- 2
#' g2 <- 2
#' chisq_g(x, g1, g2)
#'
#' x <- c(10, 25, 10, 0, 0)
#' g1 <- 1
#' g2 <- 1
#' chisq_g(x, g1, g2)
#'
#' @export
chisq_g <- function(x, g1, g2){
  TOL <- sqrt(.Machine$double.eps)
  ploidy <- length(x) - 1
  nd <- n_dr_params(ploidy = ploidy)
  gf <- gf_freq(
    p1_g = g1,
    p1_ploidy = ploidy,
    p1_gamma = NULL,
    p1_alpha = rep(0, times = nd),
    p1_beta = NULL,
    p1_type = "polysomic",
    p1_add_dr = FALSE,
    p2_g = g2,
    p2_ploidy = ploidy,
    p2_gamma = NULL,
    p2_alpha = rep(0, times = nd),
    p2_beta = NULL,
    p2_type = "polysomic",
    p2_add_dr = FALSE,
    pi = 0,
    nudge = 0)
  which_zero <- gf < TOL
  gf[which_zero] <- 0
  gf <- gf / sum(gf)


  if (sum(x[which_zero]) > 0.5) { ## if any incompatibility, p-value is 0
    ret <- list(statistic = Inf,
                p_value = 0,
                df = NA_real_)
    return(ret)
  } else {
    x <- x[!which_zero]
    gf <- gf[!which_zero]
  }

  if (length(x) == 1) {
    ret <- list(statistic = 0, p_value = 1, df = 0)
  } else {
    suppressWarnings( ## for small sample size
      chout <- stats::chisq.test(x = x, p = gf)
    )
    ret <- list(statistic = chout$statistic[[1]],
                p_value = chout$p.value[[1]],
                df = chout$parameter[[1]])
  }
  return(ret)

}

#' @describeIn chisq_g Alias for chisq_g, for backwards compatibility.
#'
#' @export
chisq_g4 <- chisq_g

#' @describeIn chisq_gl Alias for chisq_gl, for backwards compatibility.
#'
#' @export
chisq_gl4 <- chisq_gl
