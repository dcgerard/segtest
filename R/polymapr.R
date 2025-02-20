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
#'  \item{bestfit}{The genotype frequencies of the best fit model.}
#'  \item{frq_invalid}{The frequency of invalid genotypes.}
#'  \item{p_invalid}{The p-value of the invalid proportion.}
#' }
#'
#' @seealso \code{\link[polymapR]{checkF1}()}.
#'
#' @examples
#' g1 <- 0
#' g2 <- 1
#' x <- c(4, 16, 0, 0, 0)
#' polymapr_test(x = x, g1 = g1, g2 = g2, type = "segtest")
#' polymapr_test(x = x, g1 = g1, g2 = g2, type = "polymapR")
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
  if (is.matrix(x)) {
    dat <- "gl"
  } else {
    dat <- "known"
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
  p_invalid <- stats::pbinom(
    q = (1 - cout$checked_F1$frqInvalid_bestParentfit[[1]]) * n,
    size = n,
    prob = 1 - seg_invalidrate)

  ret <- list(
    p_value = cout$checked_F1$Pvalue_bestParentfit[[1]],
    bestfit = polymapr_model_to_prop(mod = as.character(cout$checked_F1$bestParentfit[[1]]), ploidy = ploidy),
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
  ploidy <- ncol(gl) - 1
  stopifnot(abs(rowSums(gl) - 1) < sqrt(.Machine$double.eps))
  seg_invalidrate <- 0.03

  n <- nrow(gl)

  ## need to repeat parent info to force polymapR to use it
  gl <- rbind(
    gl,
    (0:ploidy == g1) * 1,
    (0:ploidy == g1) * 1,
    (0:ploidy == g2) * 1,
    (0:ploidy == g2) * 1
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
    ploidy = ploidy,
    polysomic = TRUE,
    disomic = TRUE,
    mixed = TRUE,
    F1 = paste0("F", seq_len(n))
  )

  p_invalid <- stats::pbinom(q = (1 - cout$checked_F1$frqInvalid_bestParentfit[[1]]) * n,
                             size = n,
                             prob = 1 - seg_invalidrate)

  ret <- list(
    p_value = cout$checked_F1$Pvalue_bestParentfit[[1]],
    bestfit = polymapr_model_to_prop(mod = as.character(cout$checked_F1$bestParentfit[[1]]), ploidy = ploidy),
    frq_invalid = cout$checked_F1$frqInvalid_bestParentfit[[1]],
    p_invalid = p_invalid
  )

  return(ret)
}

#' Test for segregation distortion, approximating polymapR procedure
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
  ploidy <- length(x) - 1
  n <- sum(x)
  stopifnot(ploidy %% 2 == 0, x >= 0)
  if (!is.null(g1)) {
    stopifnot(g1 >= 0, g1 <= ploidy)
  }
  if (!is.null(g2)) {
    stopifnot(g2 >= 0, g2 <= ploidy)
  }

  if (is.null(g1)) {
    p1_pos <- 0:ploidy
  } else {
    p1_pos <- g1
  }
  if (is.null(g2)) {
    p2_pos <- 0:ploidy
  } else {
    p2_pos <- g2
  }

  TOL <- sqrt(.Machine$double.eps)
  pval <- 0
  p_invalid <- 0
  q_best <- NULL
  frq_invalid <- NULL
  for (p1_geno in p1_pos) {
    p1_list <- segtest::seg[segtest::seg$ploidy == ploidy & segtest::seg$g == p1_geno, ]$p
    for (i in seq_along(p1_list)) {
      for (p2_geno in p2_pos) {
        p2_list <- segtest::seg[segtest::seg$ploidy == ploidy & segtest::seg$g == p2_geno, ]$p
        for (j in seq_along(p2_list)) {
          fq <- convolve_2(p1 = p1_list[[i]], p2 = p2_list[[j]], nudge = 0)
          not_0 <- fq > sqrt(.Machine$double.eps)
          if (any(x[not_0] != 0)) {
            if (sum(not_0) == 1) {
              chout <- list(p.value = 1)
            } else {
              suppressWarnings( ## small sample size warning
                chout <- stats::chisq.test(x = x[not_0], p = fq[not_0])
              )
            }
            chout$frq_invalid <- sum(x[!not_0])
            chout$p_invalid <- stats::pbinom(
              q = n - chout$frq_invalid,
              size = n,
              prob = 1 - seg_invalidrate
              )

            ## weird criterion
            if (pval * p_invalid < chout$p.value * chout$p_invalid) {
              pval <- chout$p.value
              frq_invalid <- chout$frq_invalid
              q_best <- fq
              p_invalid <- chout$p_invalid
            }
          }
        }
      }
    }
  }

  ret <- list(
    p_value = pval,
    best_fit = q_best,
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
  ploidy <- ncol(gl) - 1
  stopifnot(abs(rowSums(gl) - 1) < sqrt(.Machine$double.eps))
  x <- colSums(gl)
  ## remove low counts and renormalize to sum to number of individuals
  proba_correct <- 0.05 * nrow(gl) / (ploidy + 1)
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

#' Convert polymapR's model shorthand to actual genotype frequencies
#'
#' @param mod Model of the form something like 121_2 which corrsponds to
#'    segregation ratios of 1 to 2 to 1 starting at genotype 2.
#' @param ploidy Ploidy of the offspring.
#'
#' @author David Gerard
#'
#' @examples
#' polymapr_model_to_prop(mod = "BZB_1", ploidy = 4)
#' q <- c(0, 11, 35, 11, 0)
#' q / sum(q)
#'
#' polymapr_model_to_prop(mod = "1412_0", ploidy = 4)
#' q <- c(1, 4, 1, 2, 0)
#' q / sum(q)
#'
#'
#' @noRd
polymapr_model_to_prop <- function(mod, ploidy) {
  sout <- strsplit(mod, split = "")[[1]]
  letvec <- match(x = sout, table = LETTERS)
  if (any(!is.na(letvec))) {
    sout[!is.na(letvec)] <- (10:35)[letvec[!is.na(letvec)]]
  }
  un <- which(sout == "_")
  rat <- as.numeric(sout[1:(un - 1)])
  st <- as.numeric(sout[(un + 1):length(sout)])
  fq <- rep(0, times = ploidy + 1)
  fq[(st + 1):(st + length(rat))] <- rat
  fq <- fq / sum(fq)
  return(fq)
}
