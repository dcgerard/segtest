## polymapR's version of checking for segregation distortion


#' Run segregation distortion tests as implemented in the polymapR package.
#'
#' The polymapR package tests for segregation distortion by iterating through all
#' possible forms of disomic or polysomic inheritance from either parent,
#' tests for concordance of the offspring genotypes using a chi-squared
#' test, and returns the largest p-value. It sometimes chooses a different
#' p-value based on other heuristics. They also sometimes return NA.
#' When \code{type = "segtest"}, we only look at patterns of the
#' given parent genotypes, choosing the largest p-value. When
#' \code{type = "polymapR"}, we return what they use via their heuristics.
#'
#'
#' @param x Either a vector of genotype counts, or a matrix of genotype
#'     posteriors where the rows index the individuals and the columns
#'     index the genotypes.
#' @param g1 Parent 1's genotype.
#' @param g2 Parent 2's genotype.
#' @param type Either my implementation which approximates that of
#'     polymapR (\code{"segtest"}) or the implementation
#'     through polymapR (\code{"polymapR"}). Note that
#'     polymapR needs to be installed for \code{type = "polymapR"}.
#'
#' @return A list with the following elements:
#' \describe{
#'  \item{p_value}{The p-value of the test.}
#'  \item{bestfit}{The best fit model, using the same notation as in
#'      \code{\link[polymapR]{checkF1}()}.}
#'  \item{frq_invalid}{The frequency of invalid genotypes.}
#' }
#'
#' @seealso \code{\link[polymapR]{checkF1}()}.
#'
#' @examples
#' g1 <- 0
#' g2 <- 1
#' x <- c(4, 16, 0, 0, 0)
#' polymapr_test(x = x, g1 = g1, g2 = g2, type = "segtest")
#'
#'
#'
#' @author David Gerard
#'
#' @export
polymapr_test <- function(x, g1 = NULL, g2 = NULL, type = c("segtest", "polymapR")) {
  type <- match.arg(type)
  if (!requireNamespace("polymapR", quietly = TRUE) && type == "polymapR") {
    stop(
      paste0(
        "polymapr_test:\n",
        "polymapR needs to be installed to use type = 'polymapR'.\n",
        "You can install it with install.packages('polymapR')"
        )
      )
  }

  if ((is.null(g1) | is.null(g2)) && type == "polymapR") {
    stop("unknown parent genotypes only for type = 'segtest'")
  }

  dat <- NULL
  if (length(x) == 5) {
    dat <- "known"
  } else if (ncol(x) == 5) {
    dat <- "gl"
  } else {
    stop("x needs to either have length 5 or be a matrix with 5 columns.")
  }

  if (dat == "known" && type == "segtest") {
    ret <- polymapr_approx_g(x = x, g1 = g1, g2 = g2)
  } else if (dat == "gl" && type == "segtest") {
    ret <- polymapr_approx_gl(gl = x, g1 = g1, g2 = g2)
  } else if (dat == "known" && type == "polymapR") {
    ret <- polymapR_package_g(x = x, g1 = g1, g2 = g2)
  } else if (dat == "gl" && type == "polymapR") {
    ret <- polymapr_package_gl(gl = x, g1 = g1, g2 = g2)
  } else {
    stop("dud")
  }

  return(ret)
}


#' polymapR test when genotypes are known.
#'
#' @param x genotype count vector
#' @param g1 parent 1's gentoype
#' @param g2 parent 2's genotype
#'
#' @author David Gerard
#'
#' @examples
#' x <- c(3, 22, 50, 22, 3)
#' polymapR_package_g(x = x, g1 = 2, g2 = 2)
#'
#' @noRd
polymapR_package_g <- function(x, g1, g2) {
  stopifnot(requireNamespace("polymapR", quietly = TRUE))
  seg_invalidrate <- 0.03
  ploidy <- length(x) - 1
  stopifnot(
    g1 >= 0,
    g2 >= 0,
    g1 <= ploidy,
    g2 <= ploidy
  )

  n <- sum(x)
  ## need repetitions of parent genotypes to force polymapR to keep those as info
  df <- matrix(c(g1, g1, g2, g2, gcount_to_gvec(gcount = x)), nrow = 1)
  fnames <- paste0("F", seq_len(ncol(df) - 4))
  colnames(df) <- c("P11", "P12", "P21", "P22", fnames)

  df <- rbind(df, 1) ## they forgot a drop = FALSE somewhere

  rownames(df) <- c("M1", "M2")

  cout <- polymapR::checkF1(
    input_type = "discrete",
    dosage_matrix = df,
    parent1 = c("P11", "P12"),
    parent2 = c("P21", "P22"),
    F1 = fnames,
    polysomic = TRUE,
    disomic = TRUE,
    mixed = TRUE,
    ploidy = ploidy)

  TOL <- sqrt(.Machine$double.eps) ## to get >= in pbinom
  p_invalid <- stats::pbinom(q = (1 - cout$checked_F1$frqInvalid_bestParentfit[[1]]) * n,
                             size = n,
                             prob = 1 - seg_invalidrate)

  ret <- list(
    p_value = cout$checked_F1$Pvalue_bestParentfit[[1]],
    bestfit = as.character(cout$checked_F1$bestParentfit[[1]]),
    frq_invalid = cout$checked_F1$frqInvalid_bestParentfit[[1]],
    p_invalid = p_invalid
  )

  return(ret)
}

#' @param gl posterior probability matrix. Rows index individuals, columns
#'     index genotypes
#'
#' @author David Gerard
#'
#' @examples
#' g1 <- 1
#' g2 <- 2
#' gl <- simf1gl(n = 100, g1 = g1, g2 = g2)
#' gl <- exp(gl - apply(gl, 1, log_sum_exp))
#'
#'
#' @noRd
polymapr_package_gl <- function(gl, g1, g2) {
  stopifnot(ncol(gl) == 5)
  stopifnot(abs(rowSums(gl) - 1) < sqrt(.Machine$double.eps))
  seg_invalidrate <- 0.03

  n <- nrow(gl)

  ## need to repeat parent info to force polymapR to use it
  gl <- rbind(
    gl,
    (0:4 == g1) * 1,
    (0:4 == g1) * 1,
    (0:4 == g2) * 1,
    (0:4 == g2) * 1
  )

  gp_df <- as.data.frame(gl)
  colnames(gp_df) <- paste0("P", seq_len(ncol(gl)) - 1)
  gp_df$SampleName <- c(paste0("F", seq_len(n)), "P11", "P12", "P21", "P22")
  gp_df$MarkerName <- "M1"
  gp_df$maxP <- apply(gl, 1, max)
  gp_df$maxgeno <- apply(gl, 1, which.max) - 1
  gp_df$geno <- gp_df$maxgeno

  cout <- polymapR::checkF1(
    input_type = "probabilistic",
    probgeno_df = gp_df,
    parent1 = c("P11", "P12"),
    parent2 = c("P21", "P22"),
    ploidy = 4,
    polysomic = TRUE,
    disomic = TRUE,
    mixed = TRUE,
    F1 = paste0("F", seq_len(n))
  )

  TOL <- sqrt(.Machine$double.eps) ## to get >= in pbinom
  p_invalid <- stats::pbinom(q = (1 - cout$checked_F1$frqInvalid_bestParentfit[[1]]) * n,
                             size = n,
                             prob = 1 - seg_invalidrate)

  ret <- list(
    p_value = cout$checked_F1$Pvalue_bestParentfit[[1]],
    bestfit = as.character(cout$checked_F1$bestParentfit[[1]]),
    frq_invalid = cout$checked_F1$frqInvalid_bestParentfit[[1]],
    p_invalid = p_invalid
  )

  return(ret)
}

#' test for segregation distortion, approximating polymapR procedure
#'
#' @inheritParams polymapR
#' @param seg_invalidrate If there is only one class possible, the p-value
#'     is the binomial probability of the invalid number against a true
#'     invalid rate of this, defaults to 0.03.
#'
#' @author David Gerard
#'
#' @noRd
polymapr_approx_g <- function(x, g1 = NULL, g2 = NULL, seg_invalidrate = 0.03) {
  ploidy <- 4
  n <- sum(x)
  stopifnot(length(x) == 5, x >= 0)
  if (!is.null(g1)) {
    stopifnot(g1 >= 0, g1 <= ploidy)
  }
  if (!is.null(g2)) {
    stopifnot(g2 >= 0, g2 <= ploidy)
  }

  TOL <- sqrt(.Machine$double.eps)
  pval <- 0
  p_invalid <- 0
  bi <- NULL
  frq_invalid <- NULL
  for (i in seq_len(nrow(segtypes))) {
    if (!is.null(g1) && !is.null(g2)) {
      is_p <- any((segtypes$pardosage[[i]][, 1] == g1) &
                    (segtypes$pardosage[[i]][, 2] == g2))
    } else if (!is.null(g1)) {
      is_p <- any((segtypes$pardosage[[i]][, 1] == g1))
    } else if (!is.null(g2)) {
      is_p <- any((segtypes$pardosage[[i]][, 2] == g2))
    } else {
      is_p <- TRUE
    }
    if (is_p) {
      fq <- segtypes$freq[[i]]
      not_0 <- fq > sqrt(.Machine$double.eps)
      if (any(x[not_0] != 0)) {
        if (sum(not_0) == 1) {
          chout <- list(p.value = 1)
        } else {
          suppressWarnings(
            chout <- stats::chisq.test(x = x[not_0], p = fq[not_0])
          )
        }
        chout$frq_invalid <- sum(x[!not_0])
        chout$p_invalid <- stats::pbinom(q = n - chout$frq_invalid,
                                         size = n,
                                         prob = 1 - seg_invalidrate)

        ## weird criterion
        if (pval * p_invalid < chout$p.value * chout$p_invalid) {
          pval <- chout$p.value
          bi <- i
          frq_invalid <- chout$frq_invalid
          p_invalid <- chout$p_invalid
        }
      }
    }
  }

  ret <- list(
    p_value = pval,
    best_fit = segtypes$mod[[bi]],
    frq_invalid = frq_invalid,
    p_invalid = p_invalid
  )
  return(ret)
}

#' @param gl posterior probability matrix. Rows index individuals, columns
#'     index genotypes
#'
#' @author David Gerard
#'
#' @noRd
polymapr_approx_gl <- function(gl, g1, g2, seg_invalidrate = 0.03) {
  stopifnot(ncol(gl) == 5)
  stopifnot(abs(rowSums(gl) - 1) < sqrt(.Machine$double.eps))
  x <- colSums(gl)
  ## remove low counts and renormalize to sum to number of individuals
  proba_correct <- 0.05 * nrow(gl) / (4 + 1)
  x[x < proba_correct] <- 0
  x <- (x/sum(x)) * nrow(gl)
  ## now just run it things like genotypes are known
  ret <- polymapr_approx_g(
    x = x,
    g1 = g1,
    g2 = g2,
    seg_invalidrate = seg_invalidrate)
  return(ret)
}
