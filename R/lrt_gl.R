#' Likelihood ratio test when only have genotype likelihoods
#'
#' @param gl genotype log-likelihoods of offspring. Rows index individuals
#'     and columns index genotypes.
#' @param p1_gl genotype log-likelihoods of parent 1.
#' @param p2_gl genotype log-likelihoods of parent 2.
#' @param drbound maximum double reduction rate
#' @param ntry number of restarts of optimizatioin
#' @param dr Should we allow for double reduction (\code{TRUE})
#'     or not (\code{FALSE})?
#' @param pp Should we allow for partial preferential pairing
#'     (\code{TRUE}) or not (\code{FALSE})?
#' @param alpha If \code{dr = FALSE}, this is the known rate of double
#'     reduction.
#' @param xi1 If \code{pp = FALSE}, this is the known preferential pairing
#'     parameter of parent 1.
#' @param xi2 If \code{pp = FALSE}, this is the known preferential pairing
#'     parameter of parent 2.
#'
#' @author David Gerard
#'
#' @noRd
lrt_gl4 <- function(
    gl,
    p1_gl = rep(-log(5), 5),
    p2_gl = rep(-log(5), 5),
    drbound = 1/6,
    ntry = 5,
    dr = TRUE,
    pp = TRUE,
    alpha = 0,
    xi1 = 1/3,
    xi2 = 1/3) {
  stopifnot(length(p1_gl) == 5,
            length(p2_gl) == 5,
            ncol(gl) == 5,
            length(drbound) == 1,
            drbound > 1e-6,
            drbound < 1,
            ntry >= 1,
            length(ntry) == 1,
            length(dr) == 1,
            is.logical(dr),
            length(pp) == 1,
            is.logical(pp))

  ## normalize to sum to one
  p1_gl <- p1_gl - log_sum_exp(p1_gl)
  p2_gl <- p2_gl - log_sum_exp(p2_gl)

  ## iterate over all possible parent genos, choose minimum LRT statistic
  pdf <- expand.grid(p1 = 0:4, p2 = 0:4)
  bout <- NULL
  ll <- Inf
  p1 <- NULL
  p2 <- NULL
  for (i in seq_len(nrow(pdf))) {
    if (dr && pp) {
      lout <- lrt_dr_pp_glpknown4(
        gl = gl,
        g1 = pdf$p1[[i]],
        g2 = pdf$p2[[i]],
        drbound = drbound,
        ntry = ntry)
    } else if (dr && !pp) {
      lout <- lrt_dr_npp_glpknown4(
        gl = gl,
        g1 = pdf$p1[[i]],
        g2 = pdf$p2[[i]],
        drbound = drbound,
        ntry = ntry,
        xi1 = xi1,
        xi2 = xi2)
    } else if (!dr && pp) {
      lout <- lrt_ndr_pp_glpknown4(
        gl = gl,
        g1 = pdf$p1[[i]],
        g2 = pdf$p2[[i]],
        alpha = alpha)
    } else {
      lout <- lrt_ndr_npp_glpknown4(
        gl = gl,
        g1 = pdf$p1[[i]],
        g2 = pdf$p2[[i]],
        alpha = alpha,
        xi1 = xi1,
        xi2 = xi2)
    }

    llnew <- lout$statistic - 2 * (p1_gl[[pdf$p1[[i]] + 1]] + p2_gl[[pdf$p2[[i]] + 1]])
    if (llnew < ll) {
      ll <- llnew
      bout <- c(
        lout,
        p1 = pdf$p1[[i]],
        p2 = pdf$p2[[i]])
    }
  }
  return(bout)
}
