## Test using polymapR strategy of testing both compatible
## and number of incompatible.

#' Returns a vector of length 5 which is TRUE if that genotype is possible
#' and FALSE if it is impossible
#'
#' @noRd
which_possible <- function(g1, g2, dr = TRUE) {
  if (g1 == 0 && g2 == 0) {
    return(c(TRUE, FALSE, FALSE, FALSE, FALSE))
  } else if (g1 == 4 & g2 == 4) {
    return(c(FALSE, FALSE, FALSE, FALSE, TRUE))
  } else if (g1 == 0 && g2 == 4 || g1 == 4 && g2 == 0) {
    return(c(FALSE, FALSE, TRUE, FALSE, FALSE))
  } else if (dr) {
    if (g1 == 0 && g2 %in% c(1, 2, 3) || g1 %in% c(1, 2, 3) && g2 == 0) {
      return(c(TRUE, TRUE, TRUE, FALSE, FALSE))
    } else if (g1 == 4 && g2 %in% c(1, 2, 3) || g1 %in% c(1, 2, 3) && g2 == 4) {
      return(c(FALSE, FALSE, TRUE, TRUE, TRUE))
    }
  } else if (!dr) {
    if (g1 == 0 && g2 == 1 || g1 == 1 && g2 == 0) {
      return(c(TRUE, TRUE, FALSE, FALSE, FALSE))
    } else if (g1 == 0 && g2 == 2 || g1 == 2 && g2 == 0 || g1 == 1 && g2 == 1) {
      return(c(TRUE, TRUE, TRUE, FALSE, FALSE))
    } else if (g1 == 0 && g2 == 3 || g1 == 3 && g2 == 0) {
      return(c(FALSE, TRUE, TRUE, FALSE, FALSE))
    } else if (g1 == 1 && g2 == 4 || g1 == 4 && g2 == 1) {
      return(c(FALSE, FALSE, TRUE, TRUE, FALSE))
    } else if (g1 == 2 && g2 == 4 || g1 == 4 && g2 == 2 || g1 == 3 && g2 == 3) {
      return(c(FALSE, FALSE, TRUE, TRUE, TRUE))
    } else if (g1 == 3 && g2 == 4 || g1 == 4 && g2 == 3) {
      return(c(FALSE, FALSE, FALSE, TRUE, TRUE))
    } else if (g1 == 1 && g2 == 2 || g1 == 2 && g2 == 1) {
      return(c(TRUE, TRUE, TRUE, TRUE, FALSE))
    } else if (g1 == 1 && g2 == 3 || g1 == 3 && g2 == 1) {
      return(c(FALSE, TRUE, TRUE, TRUE, FALSE))
    } else if (g1 == 2 && g2 == 3 || g1 == 3 && g2 == 2) {
      return(c(FALSE, TRUE, TRUE, TRUE, TRUE))
    }
  }

  return(c(TRUE, TRUE, TRUE, TRUE, TRUE))
}

#' Jointly tests for segregation distortion and number of incompatible genotypes
#'
#' This is experimental. I haven't tested it out in lots of scenarios yet.
#'
#' Here, we test if the compatible genotypes are consistent with F1 populations
#' and separately test that the number of incompatible genotypes isn't too
#' large (less than 3 percent by default). This is the strategy the
#' polymapR software uses. But we use a Bonferroni correction to combine
#' these tests (minimum of two times the p-values), while they just multiply
#' the p-values together. So our approach accounts for double reduction and
#' preferential pairing, while also controlling the family-wise error rate.
#'
#' @inheritParams lrt_men_g4
#' @param pbad The upper bound on the number of bad genotypes
#'
#' @examples
#' # Run a test where genotypes 0, 1, and 2 are possible
#' x <- c(10, 10, 4, 0, 5)
#' otest_g(x = x, g1 = 1, g2 = 0)
#'
#' # polymapR's multiplication and the Bonferroni differ
#' df <- expand.grid(p1 = seq(0, 1, length.out = 20), p2 = seq(0, 1, length.out = 20))
#' df$polymapr <- NA
#' df$bonferroni <- NA
#' for (i in seq_len(nrow(df))) {
#'   df$polymapr[[i]] <- df$p1[[i]] * df$p2[[i]]
#'   df$bonferroni[[i]] <- 2 * min(c(df$p1[[i]], df$p2[[i]], 0.5))
#' }
#' graphics::plot(df$polymapr, df$bonferroni)
#'
#'
#' @return A list with the following elements
#' \describe{
#'   \item{\code{statistic}}{The log-likelihood ratio test statistic.}
#'   \item{\code{df}}{The degrees of freedom.}
#'   \item{\code{p_value}}{The Bonferroni corrected p-value.}
#'   \item{\code{p_lrt}}{The p-value of the LRT.}
#'   \item{\code{p_binom}}{The p-value of the one-sided binomial test.}
#'   \item{\code{alpha}}{The estimated double reduction rate.}
#'   \item{\code{xi1}}{The estimated preferential pairing parameter of parent 1.}
#'   \item{\code{xi2}}{The estimated preferential pairing parameter of parent 2.}
#' }
#'
#' @author David Gerard
#'
#' @export
otest_g <- function(
    x,
    g1,
    g2,
    pbad = 0.03,
    drbound = 1/6,
    pp = TRUE,
    dr = TRUE,
    alpha = 0,
    xi1 = 1/3,
    xi2 = 1/3) {
  pos <- which_possible(g1 = g1, g2 = g2, dr = dr)

  y <- x
  y[!pos] <- 0
  lout <- lrt_men_g4(
    x = y,
    g1 = g1,
    g2 = g2,
    drbound = drbound,
    pp = pp,
    dr = dr,
    alpha = alpha,
    xi1 = xi1,
    xi2 = xi2)

  bout <- stats::binom.test(
    x = sum(x[!pos]),
    n = sum(x),
    p = pbad,
    alternative = "greater")

  lout$p_binom <- bout$p.value

  ## Bonferroni correction
  lout$p_lrt <- lout$p_value
  lout$p_value <- 2 * min(c(lout$p_binom, lout$p_lrt, 0.5))

  return(lout)
}
