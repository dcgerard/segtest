## Tolerance global variables
pkg_env <- new.env()
pkg_env$TOL_small <- sqrt(.Machine$double.eps) ## for everything but df calculations.
pkg_env$TOL_big <- .Machine$double.eps^(1/4) ## only used for df calculations

logit <- function (x) {
    log(x) - log(1 - x)
}

expit <- function (x) {
    1/(1 + exp(-x))
}

all_multinom <- function (n, k) {
  if (k == 0) {
    return(NULL)
  }
  else if (k == 1) {
    return(n)
  }
  mat_list <- list()
  for (i in 0:n) {
    mat_list[[i + 1]] <- cbind(i, all_multinom(n - i, k - 1))
  }
  mat <- do.call(what = rbind, args = mat_list)
  colnames(mat) <- NULL
  mat
}

#' Upper bounds on double reduction rates.
#'
#' Provides the upper bounds on the double reduction rates based on the
#' formulas in Huang et al. (2019). There are two upper bounds provided.
#' The upper bound from complete equational separation is higher than
#' the upper bound from the pure random chromatid segregation.
#'
#' @param ploidy The ploidy
#' @param model Either complete equational segregation (\code{"ces"}) (Mather, 1935)
#'    or pure random chromatid segregation \code{"prcs"})  (Haldane, 1930). See also
#'    Huang et al. (2019).
#'
#' @return A vector of length \code{floor(ploidy / 4)} containing the upper bounds
#'    on the rates of double reduction. The \code{i}the element is the
#'    upper bound on the probability that there are \code{i} pairs of
#'    identical by double reduction alleles in a gamete.
#'
#' @references
#' \itemize{
#'   \item{Haldane, J. B. S. (1930). Theoretical genetics of autopolyploids. \emph{Journal of genetics}, 22, 359-372.}
#'   \item{Huang, K., Wang, T., Dunn, D. W., Zhang, P., Cao, X., Liu, R., & Li, B. (2019). Genotypic frequencies at equilibrium for polysomic inheritance under double-reduction. \emph{G3: Genes, Genomes, Genetics}, 9(5), 1693-1706.}
#'   \item{Mather, K. (1935). Reductional and equational separation of the chromosomes in bivalents and multivalents. \emph{Journal of genetics}, 30, 53-78.}
#' }
#'
#' @author David Gerard
#'
#' @export
drbounds <- function (ploidy, model = c("ces", "prcs")) {
  model <- match.arg(model)
  stopifnot(ploidy > 2)
  stopifnot(ploidy%%2 == 0)
  ibdr <- floor(ploidy/4)
  if (model == "ces") {
    alpha <- rep(NA_real_, length = ibdr)
    for (i in seq_along(alpha)) {
        jvec <- i:ibdr
        alpha[[i]] <- exp(
          log_sum_exp(
            (ploidy/2 - 3 * jvec) * log(2) + lchoose(jvec, i) +
              lchoose(ploidy/2, jvec) +
              lchoose(ploidy/2 - jvec, ploidy/2 - 2 * jvec) -
              lchoose(ploidy, ploidy/2)
            )
          )
    }
  } else if (model == "prcs") {
    ivec <- seq_len(ibdr)
    alpha <- exp(
      lchoose(ploidy, ploidy / 2 - ivec) + lchoose(ploidy / 2 - ivec,  ivec) +
        (ploidy / 2 -  2 * ivec) * log(2) - lchoose(2 * ploidy, ploidy / 2)
    )
  }
  return(alpha)
}

#' Bounds on the distortion at simplex loci caused by double reduction.
#'
#' Returns the upper bound on the probability of a gamete with a genotype of
#' 2 when the parent has a genotype of 1. This is based on two models.
#' The upper bound from complete equational separation is higher than
#' the upper bound from the pure random chromatid segregation. See
#' Huang et al (2019) for a description of these models.
#'
#'
#' @inheritParams drbounds
#'
#' @author David Gerard
#'
#' @references
#' \itemize{
#'   \item{Haldane, J. B. S. (1930). Theoretical genetics of autopolyploids. \emph{Journal of genetics}, 22, 359-372.}
#'   \item{Huang, K., Wang, T., Dunn, D. W., Zhang, P., Cao, X., Liu, R., & Li, B. (2019). Genotypic frequencies at equilibrium for polysomic inheritance under double-reduction. \emph{G3: Genes, Genomes, Genetics}, 9(5), 1693-1706.}
#'   \item{Mather, K. (1935). Reductional and equational separation of the chromosomes in bivalents and multivalents. \emph{Journal of genetics}, 30, 53-78.}
#' }
#'
#' @export
beta_bounds <- function(ploidy, model = c("ces", "prcs")) {
  model <- match.arg(model)
  dvec <- drbounds(ploidy = ploidy, model = model)
  return(sum(seq_along(dvec) * dvec) / ploidy)
}

#' Simplex to real function
#'
#' @param q A point on the simplex
#'
#' @return A point on the reals, of length one less than q.
#'
#' @references
#' \itemize {
#'  \item{Betancourt, M. J. (2010). Cruising the simplex: Hamiltonian Monte Carlo and the Dirichlet distribution. arXiv preprint arXiv:1010.3436.}
#'  \item{\url{https://mc-stan.org/docs/2_19/reference-manual/simplex-transform-section.html}}
#' }
#'
#' @author David Gerard
#'
#' @examples
#' simplex_to_real(c(0.5, 0.3, 0.2))
#' real_to_simplex(simplex_to_real(c(0.5, 0.3, 0.2)))
#'
#' @noRd
simplex_to_real <- function (q) {
    K <- length(q)
    stopifnot(K > 1)
    z <- q/(1 - c(0, cumsum(q)[-K]))
    y <- logit(z[-K]) - log(1/(K - 1:(K - 1)))
    return(y)
}

#' Real to simplex function
#'
#' @param y A point on the reals
#'
#' @return A point on the simplex, of length one greater than y.
#'
#' @references
#' \itemize {
#'  \item{Betancourt, M. J. (2010). Cruising the simplex: Hamiltonian Monte Carlo and the Dirichlet distribution. arXiv preprint arXiv:1010.3436.}
#'  \item{\url{https://mc-stan.org/docs/2_19/reference-manual/simplex-transform-section.html}}
#' }
#'
#' @author David Gerard
#'
#' @examples
#' real_to_simplex(c(1, 3))
#' simplex_to_real(real_to_simplex(c(1, 3)))
#'
#' @noRd
real_to_simplex <- function (y) {
    stopifnot(length(y) > 0)
    K <- length(y) + 1
    z <- c(expit(y + log(1/(K - 1:(K - 1)))), 1)
    x <- rep(NA_real_, length.out = K)
    cumval <- 0
    for (i in seq_len(K)) {
        x[[i]] <- (1 - cumval) * z[[i]]
        cumval <- cumval + x[[i]]
    }
    return(x)
}

#' Gamete frequencies from meiosis under just double reduction
#'
#' Returns the probability of a gamete dosage given the parent dosage
#' (\code{g}), the parent ploidy (\code{ploidy}), and the double reduction
#' parameter(s) (\code{alpha}). This is for biallelic loci.
#'
#' @param alpha A numeric vector containing the double reduction parameter(s).
#'     This should be a vector of length \code{floor(ploidy/4)} where \code{alpha[i]}
#'     is the probability of exactly \code{i} pairs of IBDR alleles
#'     being in the gamete. Note that \code{sum(alpha)} should be less than
#'     1, as \code{1 - sum(alpha)} is the probability of no double reduction.
#' @param g The dosage of the parent. Should be an integer between \code{0}
#'     and \code{ploidy}.
#' @param ploidy The ploidy of the species. This should be an even positive
#'     integer.
#' @param log_p A logical. Should we return the log-probability (\code{TRUE})
#'     or not (\code{FALSE})? Defaults to \code{FALSE}.
#'
#' @return A vector of length \code{ploidy/2 + 1}, containing the (log)
#'     probabilities of a gamete carrying a dosage of 0,1,...,ploidy/2+1 from a
#'     parent of dosage \code{G} who has ploidy \code{ploidy} and a
#'     double reduction rate \code{alpha}.
#'
#' @author David Gerard
#'
#' @references
#' \itemize{
#' \item{Fisher, R. A., & Mather, K. (1943). The inheritance of style length in Lythrum salicaria. Annals of Eugenics, 12(1), 1-23.}
#' \item{Huang, K., Wang, T., Dunn, D. W., Zhang, P., Cao, X., Liu, R., & Li, B. (2019). Genotypic frequencies at equilibrium for polysomic inheritance under double-reduction. G3: Genes, Genomes, Genetics, 9(5), 1693-1706.}
#' }
#'
#' @noRd
#'
#' @examples
#' gamfreq_dr(alpha = 1/6, g = 2, ploidy = 4)
gamfreq_dr <- function (alpha, g, ploidy, log_p = FALSE) {
    x <- 0:(ploidy/2)
    stopifnot(length(ploidy) == 1L, length(g) == 1L, length(log_p) == 1L)
    stopifnot(is.logical(log_p))
    stopifnot(ploidy%%2 == 0)
    stopifnot(0 <= g, g <= ploidy)
    stopifnot(0 <= x, x <= ploidy/2)
    ibdr <- floor(ploidy/4)
    stopifnot(length(alpha) == ibdr)
    if (ploidy == 2) {
        retvec <- stats::dhyper(
          x = x,
          m = g,
          n = ploidy - g,
          k = ploidy / 2,
          log = log_p)
        return(retvec)
    }
    if (all(abs(alpha) < pkg_env$TOL_small)) {
        retvec <- stats::dhyper(
          x = x,
          m = g,
          n = ploidy - g,
          k = ploidy / 2,
          log = log_p)
        return(retvec)
    }
    ijmat <- cbind(utils::combn(x = 0:ibdr, m = 2), matrix(rep(0:ibdr, each = 2), nrow = 2))
    jvec <- ijmat[1, ]
    ivec <- ijmat[2, ]
    alpha <- c(1 - sum(alpha), alpha)
    alphavec <- alpha[ivec + 1]
    retvec <- rep(NA_real_, length = length(x))
    for (k in seq_along(x)) {
        retvec[[k]] <- sum(exp(lchoose(g, jvec) + lchoose(g -
            jvec, x[[k]] - 2 * jvec) + lchoose(ploidy - g, ivec -
            jvec) + lchoose(ploidy - g - (ivec - jvec), ploidy/2 -
            x[[k]] - 2 * (ivec - jvec)) - lchoose(ploidy, ivec) -
            lchoose(ploidy - ivec, ploidy/2 - 2 * ivec)) * alphavec)
    }
    retvec[is.nan(retvec)] <- -Inf
    if (log_p) {
        retvec <- log(retvec)
    }
    return(retvec)
}

#' Genotype frequencies of an F1 population under just double reduction
#'
#' Returns the distribution of an offspring dosages given
#' parental dosages (\code{g1} and \code{g2}), the ploidy of the
#' species (\code{ploidy}), and the double reduction parameter(s)
#' (\code{alpha} and \code{alpha2}).
#'
#' @inheritParams gamfreq_dr
#' @param alpha2 The double reduction rate(s) for the second parent. By default,
#'     this is the same as the first parent.
#' @param g1 The dosage of parent 1. Should be an integer between \code{0}
#'     and \code{ploidy}.
#' @param g2 The dosage of parent 2. Should be an integer between \code{0}
#'     and \code{ploidy}.
#'
#' @return A vector of probabilities. The \code{i}th element is the
#'     probability that the offspring will have dosage \code{i-1}.
#'
#' @author David Gerard
#'
#' @noRd
#'
#' @examples
#' f1_gf_dr(alpha = c(0.5, 0.1), g1 = 4, g2 = 5, ploidy = 8)
#'
f1_gf_dr <- function (alpha, alpha2 = alpha, g1, g2, ploidy) {
    stopifnot(length(g1) == 1L, length(g2) == 1L, length(ploidy) ==
        1L)
    stopifnot(ploidy%%2 == 0)
    stopifnot(0 <= g1, g1 <= ploidy)
    stopifnot(0 <= g2, g2 <= ploidy)
    ibdr <- floor(ploidy / 4)
    stopifnot(length(alpha) == ibdr)
    stopifnot(length(alpha2) == ibdr)
    p1gamprob <- gamfreq_dr(
      alpha = alpha,
      g = g1,
      ploidy = ploidy,
      log_p = FALSE)
    p2gamprob <- gamfreq_dr(
      alpha = alpha2,
      g = g2,
      ploidy = ploidy,
      log_p = FALSE)
    zygdist <- stats::convolve(p1gamprob, rev(p2gamprob), type = "open")
    zygdist[zygdist < 0] <- 0
    zygdist <- zygdist/sum(zygdist)
    return(zygdist)
}


#' Returns all gamete segregation patterns under preferential pairing
#'
#' @param ploidy The ploidy
#'
#' @return A data frame with all of the gamete frequencies possible. Columsn include
#'   \describe{
#'     \item{g}{parent genotype.}
#'     \item{m}{pairing configuration.}
#'     \item{p}{gamete frequencies.}
#'     \item{ploidy}{The ploidy}
#'   }
#'
#' @author David Gerard
#'
#' @noRd
pp_seg_pats <- function(ploidy) {
  stopifnot(ploidy%%2 == 0, ploidy >= 2)
  mmat <- all_multinom(n = ploidy / 2, k = 3)
  df <- data.frame(
    ploidy = ploidy,
    g = mmat[, 2] + 2 * mmat[, 3],
    m = I(lapply(seq_len(nrow(mmat)), function(i) mmat[i,])))
  df <- df[order(df$g), ]
  df$p <- vector(mode = "list", length = nrow(df))
  rownames(df) <- NULL
  for (i in seq_len(nrow(df))) {
    df$p[[i]] <- stats::dbinom(x = 0:(ploidy/2) - df$m[[i]][[3]], size = df$m[[i]][[2]], prob = 0.5)
  }
  return(df)
}

#' Number of mixture components
#'
#' The number of disomic inheritance patterns for a given ploidy and a given
#' parental dosage. See also \code{\link{seg}} for the list of all possible
#' disomic inheritance patterns for even ploidies up to 20.
#'
#' @param g parent genotype
#' @param ploidy parent ploidy
#'
#' @return The number of mixture components.
#'
#' @export
n_pp_mix <- function(g, ploidy) {
  ifelse(g >= ploidy / 2, ploidy / 2 - ceiling(g / 2) + 1, floor(g / 2) + 1)
}

#' Gamete frequencies from mixture of pairing configurations.
#'
#' @param gamma The mixing probabilities for the pairing configurations.
#' @param g The parent gentoype.
#' @param ploidy The ploidy of the parent.
#'
#' @author David Gerard
#'
#' @return A vector of gamete frequencies.
#'
#' @noRd
gamfreq_pp <- function(gamma, g, ploidy) {
  TOL <- pkg_env$TOL_small
  stopifnot(
    length(gamma) == n_pp_mix(g = g, ploidy = ploidy),
    abs(sum(gamma) - 1) < TOL,
    gamma >= 0)
  plist <- segtest::seg[segtest::seg$ploidy == ploidy & (segtest::seg$mode == "disomic" | segtest::seg$mode == "both") & segtest::seg$g == g, ]$p
  pvec <- rep(0, length.out = ploidy / 2 + 1)
  for (i in seq_along(plist)) {
    pvec <- pvec + plist[[i]] * gamma[[i]]
  }
  return(pvec)
}

#' Gamete frequencies from a mixture of pairing configurations and polysomic
#'
#' This is the same as \code{\link{gamfreq_pp}()} but adds a mixture
#' component of polysomic inheritance at the maximum double reduction rate.
#' The last mixture component is the one for polysomic inheritance, with
#' the effects of double reduction.
#'
#' @inheritParams gamfreq_pp
#' @param alpha The double reduction rate of the polysomic component.
#'
#' @author David Gerard
#'
#' @examples
#' gamfreq_seg(gamma = c(0.1, 0.3, 0.6), g = 2, ploidy = 6, alpha = 0.3)
#'
#'
#' @noRd
gamfreq_seg <- function(gamma, g, ploidy, alpha) {
  TOL <- pkg_env$TOL_small
  stopifnot(
    length(gamma) == n_pp_mix(g = g, ploidy = ploidy) + 1,
    abs(sum(gamma) - 1) < TOL,
    gamma >= 0)
  stopifnot(alpha >= 0, length(alpha) == floor(ploidy / 4))
  plist <- segtest::seg[segtest::seg$ploidy == ploidy & (segtest::seg$mode == "disomic" | segtest::seg$mode == "both") & segtest::seg$g == g, ]$p
  pvec <- rep(0, length.out = ploidy / 2 + 1)
  for (i in seq_along(plist)) {
    pvec <- pvec + plist[[i]] * gamma[[i]]
  }
  pvec <- pvec + gamfreq_dr(alpha = alpha, g = g, ploidy = ploidy, log_p = FALSE) * gamma[[length(gamma)]]
  return(pvec)
}

#' Gamete frequencies under a generalized model
#'
#' Returns the gamete frequencies for autopolyploids, allopolyploids,
#' and segmental allopolyploids, accounting for the effects of double
#' reduction and partial preferential pairing.
#'
#' @section Models:
#'
#' If \code{type = "polysomic"}, then the gamete frequencies correspond
#' to those of Huang et al (2019). Those formulas are for general multiallelic
#' loci, so see also Appendix G of Gerard (2022) for special case of
#' biallelic loci. The relevant parameter is \code{alpha}, a vector of
#' length \code{floor(ploidy / 4)}, where \code{alpha[[i]]} is the
#' probability that there are \code{i} pairs of double reduced alleles
#' in a gamete. The theoretical upper bound on \code{alpha} is given in
#' \code{\link{drbounds}()}.
#'
#' If \code{type = "mix"} and \code{add_dr = FALSE}, then the gamete
#' frequencies correspond to the pairing configuration model of
#' Gerard et al (2018). This model states that the gamete frequencies are
#' a convex combination of the disomic inheritance frequencies. The weights
#' of this convex combination are provided in the \code{gamma} parameter. The
#' total number of disomic segregation patterns is given by
#' \code{\link{n_pp_mix}()}. The order of these segregation patterns used is
#' the order in \code{\link{seg}}.
#'
#' The model for \code{type = "mix"} and \code{add_dr = TRUE} is the same
#' as for \code{type = "mix"} and \code{add_dr = FALSE} \emph{except} at
#' parental simplex loci. At such loci, there are no effects of preferential
#' pairing, and so the option \code{add_dr = TRUE} allows for the effects
#' of double reduction at simplex loci. The relevant parameter here is
#' \code{beta}. The first three gamete frequencies at simplex loci are
#' \code{c(0.5 + beta, 0.5 - 2 * beta, beta)}, and the rest are 0. The
#' upper bound on beta for two different models are given by
#' \code{\link{beta_bounds}()}.
#'
#' @param g Parent genotype.
#' @param ploidy Parent ploidy. Should be even, and between 2 and 20 (inclusive).
#'    Let me know if you need the ploidy to be higher. I can update
#'    the package really easily.
#' @param gamma The mixture proportions for the pairing configurations.
#'    The proportions are in the same order the configurations in
#'    \code{\link{seg}}. See Gerard et al (2018) for details on
#'    pairing configurations.
#' @param alpha The double reduction rate(s) (if using). Defaults to 0's.
#' @param beta The double reduction adjustment for simplex markers if
#'    \code{type = "mix"} and \code{add_dr = TRUE}. Assumed to be 0 by default.
#' @param type Either \code{"mix"}, meaning a mixture model of
#'    pairing configurations, or \code{"polysomic"} for polysomic inheritance.
#' @param add_dr A logical. If \code{type = "polysomic"}, then the double
#'    double reduction rate (\code{"alpha"}) will be used no matter the
#'    value of \code{add_dr}, so set \code{alpha = 0} if you don't want it.
#'    But if \code{type = "mix"} then we will incorporate double reduction
#'    only in simplex markers (where it matters the most, and where
#'    preferential pairing does not operate).
#'
#' @references
#' \itemize{
#'   \item{Gerard, D. (2023). Double reduction estimation and equilibrium tests in natural autopolyploid populations. \emph{Biometrics}, 79(3), 2143-2156.}
#'   \item{Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018). Genotyping polyploids from messy sequencing data. \emph{Genetics}, 210(3), 789-807.}
#'   \item{Huang, K., Wang, T., Dunn, D. W., Zhang, P., Cao, X., Liu, R., & Li, B. (2019). Genotypic frequencies at equilibrium for polysomic inheritance under double-reduction. \emph{G3: Genes, Genomes, Genetics}, 9(5), 1693-1706.}
#' }
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' ## Various duplex models
#' gamfreq(g = 2, ploidy = 4, gamma = c(0, 1), type = "mix")
#' gamfreq(g = 2, ploidy = 4, gamma = c(1, 0), type = "mix")
#' gamfreq(g = 2, ploidy = 4, gamma = c(0.5, 0.5), type = "mix")
#' gamfreq(g = 2, ploidy = 4, alpha = 0, type = "polysomic")
#' gamfreq(g = 2, ploidy = 4, alpha = 1/6, type = "polysomic")
#'
#' ## Various simplex models
#' gamfreq(g = 1, ploidy = 4, beta = 1/24, gamma = 1, type = "mix", add_dr = TRUE)
#' gamfreq(g = 1, ploidy = 4, alpha = 1/6, type = "polysomic")
#' gamfreq(g = 1, ploidy = 4, gamma = 1, type = "mix", add_dr = FALSE)
#' gamfreq(g = 1, ploidy = 4, alpha = 0, type = "polysomic")
gamfreq <- function(
    g,
    ploidy,
    gamma = NULL,
    alpha = NULL,
    beta = NULL,
    type = c("mix", "polysomic"),
    add_dr = TRUE) {

  if (g == 0) {
    pvec <- c(1, rep(x = 0, times = ploidy / 2))
    return(pvec)
  } else if (g == ploidy) {
    pvec <- c(rep(x = 0, times = ploidy / 2), 1)
    return(pvec)
  }

  ## Check input
  TOL <- pkg_env$TOL_small
  stopifnot(ploidy >= 2, ploidy %% 2 == 0, ploidy <= 20)
  stopifnot(g >= 0, g <= ploidy)
  if (!is.null(gamma)) {
    stopifnot(
      length(gamma) == n_pp_mix(g = g, ploidy = ploidy),
      abs(sum(gamma) - 1) < TOL,
      gamma >= 0
    )
  }
  if (!is.null(alpha)) {
    stopifnot(
      alpha >= 0,
      sum(alpha) <= 1,
      length(alpha) == floor(ploidy / 4)
    )
  }
  if (!is.null(beta)) {
    stopifnot(
      beta >= 0,
      beta <= 0.25,
      length(beta) == 1
    )
  }
  stopifnot(is.logical(add_dr), length(add_dr) == 1)
  type <- match.arg(type)
  if (type == "mix") {
    stopifnot(!is.null(gamma))
    if (add_dr) {
      if (is.null(beta)) {
        beta <- 0
      }
    }
  } else if (type == "polysomic") {
    if (is.null(alpha)) {
      alpha <- rep(0, times = floor(ploidy / 4))
    }
  }

  ## calculate gamete frequencies
  if (type == "polysomic") {
    pvec <- gamfreq_dr(alpha = alpha, g = g, ploidy = ploidy, log_p = FALSE)
  } else if (type == "mix" && add_dr && (g == 1 | g == ploidy - 1)) {
    if (g == 1) {
      pvec <- c(
        0.5 + beta,
        0.5 - 2 * beta,
        beta,
        rep(0, times = ploidy / 2 - 2)
      )
    } else if (g == ploidy - 1) {
      pvec <- c(
        rep(0, times = ploidy / 2 - 2),
        beta,
        0.5 - 2 * beta,
        0.5 + beta
      )
    }
  } else if (type == "mix") {
    pvec <- gamfreq_pp(gamma = gamma, g = g, ploidy = ploidy)
  }

  return(pvec)
}

#' Genotype frequencies of an F1 population under a generalized model.
#'
#' @param p1_g,p1_ploidy,p1_gamma,p1_alpha,p1_beta,p1_type,p1_add_dr The first parent's version of the parameters in \code{\link{gamfreq}()}.
#' @param p2_g,p2_ploidy,p2_gamma,p2_alpha,p2_beta,p2_type,p2_add_dr The second parent's version of the parameters in \code{\link{gamfreq}()}.
#' @param pi The proportion of outliers.
#'
#' @author David Gerard
#'
#' @seealso \code{\link{gamfreq}()}.
#'
#' @export
#'
#' @examples
#' q <- gf_freq(
#'   p1_g = 2,
#'   p1_ploidy = 4,
#'   p1_gamma = c(0.1, 0.9),
#'   p1_type = "mix",
#'   p2_g = 2,
#'   p2_ploidy = 6,
#'   p2_alpha = 0.1,
#'   p2_type = "polysomic",
#'   pi = 0.05)
#'
gf_freq <- function(
    p1_g,
    p1_ploidy,
    p1_gamma = NULL,
    p1_alpha = NULL,
    p1_beta = NULL,
    p1_type = c("mix", "polysomic"),
    p1_add_dr = TRUE,
    p2_g,
    p2_ploidy,
    p2_gamma = NULL,
    p2_alpha = NULL,
    p2_beta = NULL,
    p2_type = c("mix", "polysomic"),
    p2_add_dr = TRUE,
    pi = 0) {
  TOL <- pkg_env$TOL_small
  p1_type <- match.arg(p1_type)
  p2_type <- match.arg(p2_type)
  p1 <- gamfreq(
    g = p1_g,
    ploidy = p1_ploidy,
    gamma = p1_gamma,
    alpha = p1_alpha,
    beta = p1_beta,
    type = p1_type,
    add_dr = p1_add_dr)
  p2 <- gamfreq(
    g = p2_g,
    ploidy = p2_ploidy,
    gamma = p2_gamma,
    alpha = p2_alpha,
    beta = p2_beta,
    type = p2_type,
    add_dr = p2_add_dr)
  q <- stats::convolve(x = p1, y = rev(p2), type = "open")
  q[q < TOL] <- 0 ## for -1e-16
  q <- q / sum(q)
  q <- q * (1 - pi) + pi / length(q)
  return(q)
}

#' Convert parameterization from something optim() can use to something
#' gamfreq() can use.
#'
#' @param rule A list of length three. The first element gives the model
#'     of parent 1. The second element gives the model of parent 2.
#'     The third element is a logiical on whether we add an outlier or now.
#'     The first two lists have elements
#'     \describe {
#'       \item{ploidy}{The parent's ploidy.}
#'       \item{g}{The parent dosage.}
#'       \item{type}{Either "mix", "mix_dr", or "polysomic".}
#'       \item{alpha}{The fixed value of alpha (not from par}
#'       \item{beta}{The fixed value of beta (not from par)}
#'       \item{gamma}{The fixed value of gamma (not from par)}
#'     }
#'     The third list has the following elements
#'     \describe{
#'       \item{outlier}{A logical on whether to include and outlier.}
#'       \item{pi}{The fixed value of pi (not from par).}
#'     }
#' @param par A vector of parameters to be converted based on rule.
#'    \itemize{
#'    \item{
#'      If \code{type = "mix"}, then
#'      gamma is assumed to be the real-line parameterization (see
#'      \code{\link{real_to_simplex}()} and \code{\link{simplex_to_real}()})
#'      }
#'    \item{
#'      If \code{type = "mix_dr"}, then this is the same as "mix" but for simplex
#'      markers we allow for beta.
#'      Beta (if g = 1) is not real-line transformed since there are bounds on it. See
#'      \code{\link{beta_bounds}()}.
#'      }
#'    \item{
#'       If \code{type = "polysomic"}, then only alpha is specified. This is
#'       not real-line specified, since there are natural bounds on alpha.
#'       See \code{\link{drbounds}()}.
#'      }
#'    \item{
#'       All of parent 1's parameters are in the vector, then all of parent 2's
#'       parameters.
#'       }
#'    \item{
#'       If \code{outlier = TRUE}, then the outlier proportion is the last
#'       element of \code{par}
#'       }
#'    }
#'
#' @return A list of length 3. The first contains the parameters from
#'     gamfreq() for parent 1, the second contains the parameters from
#'     gamfreq of parent 2. The third contains the outlier proportion.
#'     These parameters in the parent lists are \code{ploidy}, \code{g}, \code{gamma},
#'     \code{beta}, \code{alpha}, \code{type}, and \code{add_dr}.
#'     See \code{\link{gamfreq}()} for details on these parameters. The third
#'     is a list with elements \code{outlier} (a logical on whether
#'     outliers are modeled) and \code{pi} (the outlier proportion).
#'
#' @author David Gerard
#'
#' @examples
#' rule <- list(
#'   list(ploidy = 4, g = 2, type = "mix"),
#'   list(ploidy = 8, g = 4, type = "polysomic"),
#'   list(outlier = TRUE)
#'   )
#' par <- c(-2, 0.1, 0.2, 0.03)
#' par_to_gam(par = par, rule = rule)
#'
#' rule <- list(
#'   list(ploidy = 4, g = 0, type = "mix"),
#'   list(ploidy = 8, g = 1, type = "mix"),
#'   list(outlier = TRUE)
#'   )
#' par <- c(0.1)
#' par_to_gam(par = par, rule = rule)
#'
#' @noRd
par_to_gam <- function(par, rule) {
  gam <- list()
  gam[[1]] <- list(
    ploidy = rule[[1]]$ploidy,
    g = rule[[1]]$g,
    alpha = rule[[1]]$alpha,
    beta = rule[[1]]$beta,
    gamma = rule[[1]]$gamma,
    type = NULL,
    add_dr = NULL)
  gam[[2]] <- list(
    ploidy = rule[[2]]$ploidy,
    g = rule[[2]]$g,
    alpha = rule[[2]]$alpha,
    beta = rule[[2]]$beta,
    gamma = rule[[2]]$gamma,
    type = NULL,
    add_dr = NULL)
  gam[[3]] <- list(
    outlier = rule[[3]]$outlier,
    pi = rule[[3]]$pi
  )

  ## Do parent 1 ----
  cindex <- 1 ## keeps track of where we are in par
  if (rule[[1]]$type == "polysomic") {
    gam[[1]]$type <- "polysomic"
    gam[[1]]$add_dr <- FALSE
    if (rule[[1]]$g != 0 && rule[[1]]$g != rule[[1]]$ploidy) {
      if (is.null(rule[[1]]$alpha)) {
        ndr <- floor(rule[[1]]$ploidy / 4)
        gam[[1]]$alpha <- par[cindex:(cindex + ndr - 1)]
        cindex <- cindex + ndr
      } else {
        gam[[1]]$alpha <- rule[[1]]$alpha
      }
    }
  } else if (rule[[1]]$type == "mix_dr" && (rule[[1]]$g == 1 || rule[[1]]$g == rule[[1]]$ploidy - 1)) {
    gam[[1]]$type <- "mix"
    gam[[1]]$add_dr <- TRUE
    gam[[1]]$gamma <- 1
    if (is.null(rule[[1]]$beta)) {
      gam[[1]]$beta <- par[[cindex]]
      cindex <- cindex + 1
    } else {
      gam[[1]]$beta <- rule[[1]]$beta
    }
  } else if (rule[[1]]$type == "mix" && (rule[[1]]$g == 1 || rule[[1]]$g == rule[[1]]$ploidy - 1)) {
    gam[[1]]$type <- "mix"
    gam[[1]]$gamma <- 1
    gam[[1]]$add_dr <- FALSE
  } else if ((rule[[1]]$type == "mix" || rule[[1]]$type == "mix_dr")) {
    gam[[1]]$type <- "mix"
    gam[[1]]$add_dr <- FALSE
    if (rule[[1]]$g != 0 && rule[[1]]$g != rule[[1]]$ploidy) {
      if (is.null(rule[[1]]$gamma)) {
        npp <- n_pp_mix(g = rule[[1]]$g, ploidy = rule[[1]]$ploidy)
        gam[[1]]$gamma <- real_to_simplex(par[cindex:(cindex + npp - 2)])
        cindex <- cindex + npp - 1
      } else {
        gam[[1]]$gamma <- rule[[1]]$gamma
      }
    }
  } else {
    stop("par_to_gam 1: how did you get here?")
  }

  ## Do parent 2 ----
  if (rule[[2]]$type == "polysomic") {
    gam[[2]]$type <- "polysomic"
    gam[[2]]$add_dr <- FALSE
    if (rule[[2]]$g != 0 && rule[[2]]$g != rule[[2]]$ploidy) {
      if (is.null(rule[[2]]$alpha)) {
        ndr <- floor(rule[[2]]$ploidy / 4)
        gam[[2]]$alpha <- par[cindex:(cindex + ndr - 1)]
        cindex <- cindex + ndr
      } else {
        gam[[2]]$alpha <- rule[[2]]$alpha
      }
    }
  } else if (rule[[2]]$type == "mix_dr" && (rule[[2]]$g == 1 || rule[[2]]$g == rule[[2]]$ploidy - 1)) {
    gam[[2]]$type <- "mix"
    gam[[2]]$add_dr <- TRUE
    gam[[2]]$gamma <- 1
    if (is.null(rule[[2]]$beta)) {
      gam[[2]]$beta <- par[[cindex]]
      cindex <- cindex + 1
    } else {
      gam[[2]]$beta <- rule[[2]]$beta
    }
  } else if (rule[[2]]$type == "mix" && (rule[[2]]$g == 1 || rule[[2]]$g == rule[[2]]$ploidy - 1)) {
    gam[[2]]$type <- "mix"
    gam[[2]]$gamma <- 1
    gam[[2]]$add_dr <- FALSE
  } else if ((rule[[2]]$type == "mix" || rule[[2]]$type == "mix_dr")) {
    gam[[2]]$type <- "mix"
    gam[[2]]$add_dr <- FALSE
    if (rule[[2]]$g != 0 && rule[[2]]$g != rule[[2]]$ploidy) {
      if (is.null(rule[[2]]$gamma)) {
        npp <- n_pp_mix(g = rule[[2]]$g, ploidy = rule[[2]]$ploidy)
        gam[[2]]$gamma <- real_to_simplex(par[cindex:(cindex + npp - 2)])
        cindex <- cindex + npp - 1
      } else {
        gam[[2]]$gamma <- rule[[2]]$gamma
      }
    }
  } else {
    stop("par_to_gam 2: how did you get here?")
  }

  if (rule[[3]]$outlier) {
    if (is.null(rule[[3]]$pi)) {
      gam[[3]]$pi <- par[[cindex]]
      cindex <- cindex + 1
    } else {
      gam[[3]]$pi <- rule[[3]]$pi
    }
  }

  return(gam)
}

#' Par to genotype frequencies
#'
#' @inheritParams par_to_gam
#' @param nudge how much to move above 0.
#'
#' @return The vector of genotype frequencies.
#'
#' @author David Gerard
#'
#' @examples
#' rule <- list(
#'   list(ploidy = 4, g = 2, type = "mix"),
#'   list(ploidy = 8, g = 4, type = "polysomic"),
#'   list(outlier = TRUE)
#'   )
#' par <- c(-2, 0.1, 0.2, 0.03)
#' par_to_gf(par = par, rule = rule)
#'
#' rule <- list(
#'   list(ploidy = 4, g = 0, type = "mix"),
#'   list(ploidy = 8, g = 8, type = "polysomic"),
#'   list(outlier = FALSE)
#'   )
#' par <- c()
#' par_to_gf(par = par, rule = rule)
#'
#' @noRd
par_to_gf <- function(par, rule, nudge = sqrt(.Machine$double.eps)) {
  TOL <- pkg_env$TOL_small
  gampar <- par_to_gam(par = par, rule = rule)
  p1 <- do.call(what = gamfreq, args = gampar[[1]])
  p2 <- do.call(what = gamfreq, args = gampar[[2]])
  q <- stats::convolve(x = p1, y = rev(p2), type = "open")
  q[q < TOL] <- nudge ## for -1e-16
  q <- q / sum(q)
  if (gampar[[3]]$outlier) {
    q <- q * (1 - gampar[[3]]$pi) + gampar[[3]]$pi / length(q)
  }
  return(q)
}


#' Inverse function of par_to_gam()
#'
#' @param gam A list of length 3. The first contains the parameters from
#'     gamfreq() for parent 1, the second contains the parameters from
#'     gamfreq() of parent 2. The third contains the outlier proportion.
#'     The parameters in the parent lists are
#'     \describe{
#'       \item{\code{ploidy}}{Ploidy of the parent}
#'       \item{\code{g}}{Genotype of the parent}
#'       \item{\code{gamma}}{Preferential pairing parameter of parent}
#'       \item{\code{beta}}{Double reduction adjustment at simplex markers}
#'       \item{\code{alpha}}{Double reduction rates}
#'       \item{\code{type}}{Either \code{"mix"} or \code{"polysomic"}.}
#'       \item{\code{add_dr}}{Logical. Under "mix" should we adjust for double reduction at simplex loci?}
#'     }
#'     The third list contains the following elements
#'     \describe{
#'       \item{\code{outlier}}{A logical. Should we account for outliers?}
#'       \item{\code{pi}}{The outlier proportion.}
#'     }
#' @param fix_list A list of length three, determining which elements are fixed
#'    (not part of \code{par}). They are all logicals. Each list element corresponds
#'    to a list element in \code{gam}. E.g. \code{fix_list[[1]]$alpha = TRUE}
#'    would mean that \code{gam[[1]]$alpha} is fixed (not part of \code{par}).
#'    This can be NULL if not elements are fixed.
#' @param db Bound on double reduction rate(s). Should we use the CES or the
#'    PRCS model for the upper bounds of the double reduction rate.
#' @param ob Upper bound on the outlier proportion.
#'
#'
#' @return A list of length 4, containing \code{par}, \code{rule},
#'     \code{lower}, and \code{upper}. See \code{\link{par_to_gam}()}
#'     for a description of \code{par} and \code{rule}. The \code{lower}
#'     and \code{upper} elements are the bounds on \code{par}. \code{name}
#'     tells you what each element of par corresponds to.
#'
#' @seealso \code{\link{gamfreq}()} for more details on the parent parameters.
#'
#' @author David Gerard
#'
#' @noRd
#'
#' @examples
#' gam <- list()
#' gam[[1]] <- list(
#'   ploidy = 4,
#'   g = 1,
#'   alpha = 1/6,
#'   beta = NULL,
#'   gamma = NULL,
#'   type = "polysomic",
#'   add_dr = FALSE)
#' gam[[2]] <- list(
#'   ploidy = 8,
#'   g = 4,
#'   alpha = NULL,
#'   beta = NULL,
#'   gamma = c(0.1, 0.8, 0.1),
#'   type = "mix",
#'   add_dr = FALSE)
#' gam[[3]] <- list(
#'   outlier = TRUE,
#'   pi = 0.03
#' )
#' gam_to_par(gam)
gam_to_par <- function(gam, fix_list = NULL, db = c("ces", "prcs"), ob = 0.03) {
  ret <- list()
  ret$par <- c()
  ret$rule <- list()
  ret$lower <- c()
  ret$upper <- c()
  ret$name <- c()
  db <- match.arg(db)

  ## prepopulate rule
  ret$rule[[1]] <- list(
    ploidy = gam[[1]]$ploidy,
    g = gam[[1]]$g,
    type = NULL
  )
  ret$rule[[2]] <- list(
    ploidy = gam[[2]]$ploidy,
    g = gam[[2]]$g,
    type = NULL
  )
  ret$rule[[3]] <- list(
    outlier = NULL
  )

  if (is.null(fix_list)) {
    fix_list <- list(
      list(
        alpha = FALSE,
        beta = FALSE,
        gamma = FALSE
      ),
      list(
        alpha = FALSE,
        beta = FALSE,
        gamma = FALSE
      ),
      list(
        pi = FALSE
      )
    )
  }

  ## smallest and largest values for transformed gamma parameter
  lower_val <- -1e3
  upper_val <- 1e3
  TOL <- pkg_env$TOL_small

  if (gam[[1]]$g == 0 || gam[[1]]$g == gam[[1]]$ploidy) {
    ret$rule[[1]]$type <- gam[[1]]$type
  } else if (gam[[1]]$type == "polysomic" && gam[[1]]$g != 0 && gam[[1]]$g != gam[[1]]$ploidy) {
    if (!fix_list[[1]]$alpha) {
      nd <- floor(gam[[1]]$ploidy / 4)
      stopifnot(length(gam[[1]]$alpha) == nd)
      ret$par <- c(ret$par, gam[[1]]$alpha)
      ret$lower <- c(ret$lower, rep(x = TOL, times = nd))
      ret$upper <- c(ret$upper, drbounds(ploidy = gam[[1]]$ploidy, model = db))
      ret$name <- c(ret$name, rep(x = "alpha_1", times = nd))
    } else {
      ret$rule[[1]]$alpha <- gam[[1]]$alpha
    }
    ret$rule[[1]]$type <- "polysomic"
  } else if (gam[[1]]$type == "mix" && gam[[1]]$add_dr && (gam[[1]]$g == 1 || gam[[1]]$g == gam[[1]]$ploidy - 1)) {
    if (!fix_list[[1]]$beta) {
      stopifnot(length(gam[[1]]$beta) == 1)
      ret$par <- c(ret$par, gam[[1]]$beta)
      ret$lower <- c(ret$lower, TOL)
      ret$upper <- c(ret$upper, beta_bounds(ploidy = gam[[1]]$ploidy, model = db))
      ret$name <- c(ret$name, "beta_1")
    } else {
      ret$rule[[1]]$beta <- gam[[1]]$beta
    }
    ret$rule[[1]]$type <- "mix_dr"
  } else if (gam[[1]]$type == "mix" && !gam[[1]]$add_dr && (gam[[1]]$g == 1 || gam[[1]]$g == gam[[1]]$ploidy - 1)) {
    ret$rule[[1]]$type <- "mix"
    ret$rule[[1]]$beta <- 0
  } else if (gam[[1]]$type == "mix" && gam[[1]]$g > 1 && gam[[1]]$g < gam[[1]]$ploidy - 1) {
    if (!fix_list[[1]]$gamma) {
      npp <- n_pp_mix(g = gam[[1]]$g, ploidy = gam[[1]]$ploidy)
      stopifnot(length(gam[[1]]$gamma) == npp)
      ret$par <- c(ret$par, simplex_to_real(q = gam[[1]]$gamma))
      ret$lower <- c(ret$lower, rep(x = lower_val, times = npp - 1))
      ret$upper <- c(ret$upper, rep(x = upper_val, times = npp - 1))
      ret$name <- c(ret$name, rep(x = "gamma_1", times = npp - 1))
    } else {
      ret$rule[[1]]$gamma <- gam[[1]]$gamma
    }
    ret$rule[[1]]$type <- "mix"
  }

  if (gam[[2]]$g == 0 || gam[[2]]$g == gam[[2]]$ploidy) {
    ret$rule[[2]]$type <- gam[[2]]$type
  } else if (gam[[2]]$type == "polysomic" && gam[[2]]$g != 0 && gam[[2]]$g != gam[[2]]$ploidy) {
    if (!fix_list[[2]]$alpha) {
      nd <- floor(gam[[2]]$ploidy / 4)
      stopifnot(length(gam[[2]]$alpha) == nd)
      ret$par <- c(ret$par, gam[[2]]$alpha)
      ret$lower <- c(ret$lower, rep(x = TOL, times = nd))
      ret$upper <- c(ret$upper, drbounds(ploidy = gam[[2]]$ploidy, model = db))
      ret$name <- c(ret$name, rep(x = "alpha_2", times = nd))
    } else {
      ret$rule[[2]]$alpha <- gam[[2]]$alpha
    }
    ret$rule[[2]]$type <- "polysomic"
  } else if (gam[[2]]$type == "mix" && gam[[2]]$add_dr && (gam[[2]]$g == 1 || gam[[2]]$g == gam[[2]]$ploidy - 1)) {
    if (!fix_list[[2]]$beta) {
      stopifnot(length(gam[[2]]$beta) == 1)
      ret$par <- c(ret$par, gam[[2]]$beta)
      ret$lower <- c(ret$lower, TOL)
      ret$upper <- c(ret$upper, beta_bounds(ploidy = gam[[2]]$ploidy, model = db))
      ret$name <- c(ret$name, "beta_2")
    } else {
      ret$rule[[2]]$beta <- gam[[2]]$beta
    }
    ret$rule[[2]]$type <- "mix_dr"
  } else if (gam[[2]]$type == "mix" && !gam[[2]]$add_dr && (gam[[2]]$g == 1 || gam[[2]]$g == gam[[2]]$ploidy - 1)) {
    ret$rule[[2]]$type <- "mix"
    ret$rule[[2]]$beta <- 0
  } else if (gam[[2]]$type == "mix" && gam[[2]]$g > 1 && gam[[2]]$g < gam[[2]]$ploidy - 1) {
    if (!fix_list[[2]]$gamma) {
      npp <- n_pp_mix(g = gam[[2]]$g, ploidy = gam[[2]]$ploidy)
      stopifnot(length(gam[[2]]$gamma) == npp)
      ret$par <- c(ret$par, simplex_to_real(q = gam[[2]]$gamma))
      ret$lower <- c(ret$lower, rep(x = lower_val, times = npp - 1))
      ret$upper <- c(ret$upper, rep(x = upper_val, times = npp - 1))
      ret$name <- c(ret$name, rep(x = "gamma_2", times = npp - 1))
    } else {
      ret$rule[[2]]$gamma <- gam[[2]]$gamma
    }
    ret$rule[[2]]$type <- "mix"
  }

  if (gam[[3]]$outlier) {
    if(!fix_list[[3]]$pi) {
      stopifnot(length(gam[[3]]$pi) == 1)
      ret$par <- c(ret$par, gam[[3]]$pi)
      ret$lower <- c(ret$lower, TOL)
      ret$upper <- c(ret$upper, ob)
      ret$name <- c(ret$name, "pi")
    } else {
      ret$rule[[3]]$pi <- gam[[3]]$pi
    }
    ret$rule[[3]]$outlier <- TRUE
  } else {
    ret$rule[[3]]$outlier <- FALSE
  }

  return(ret)
}


#' (log) likelihood when genotypes are known
#'
#' @inheritParams par_to_gam
#' @param x vector of counts of offspring genotypes
#' @param log_p A logical on whether we should log the likelihood or not.
#'
#' @author David Gerard
#'
#' @noRd
ll_g <- function(par, rule, x, log_p = TRUE) {
  gf <- par_to_gf(par = par, rule = rule)
  return(stats::dmultinom(x = x, prob = gf, log = log_p))
}

#' (log) likelihood when genotypes are not known
#'
#' @inheritParams par_to_gam
#' @param x A matrix of genotype log-likelihoods. The rows index the
#'     individuals and the columns index the genotypes.
#' @param log_p A logical on whether we should log the likelihood or not.
#'
#' @author David Gerard
#'
#' @noRd
ll_gl <- function(par, rule, x, log_p = TRUE) {
  gf <- par_to_gf(par = par, rule = rule)
  ll <- llike_li(B = x, lpivec = log(gf))
  return(ifelse(log_p, ll, exp(ll)))
}

#' Test for segregation distortion in a polyploid F1 population.
#'
#' @param x The data. Can be one of two forms:
#'   \itemize{
#'     \item{A vector of genotype counts. This is when offspring genotypes are known.}
#'     \item{A matrix of genotype log-likelihoods. This is when there is genotype uncertainty.
#'           The rows index the individuals and the columns index the possible genotypes.
#'           The genotype log-likelihoods should be base e (natural log).}
#'   }
#' @param p1_ploidy,p2_ploidy The ploidy of the first or second parent. Should be even.
#' @param p1,p2 One of three forms:
#'   \itemize{
#'     \item{The known genotype of the first or second parent.}
#'     \item{The vector of genotype log-likelihoods of the first or second parent. Should be base e (natural log).}
#'     \item{\code{NULL} (completely unknown)}
#'   }
#' @param model One of six forms:
#'   \describe{
#'     \item{\code{"seg"}}{Segmental allopolyploid. Allows for arbitrary levels of polysomic and disomic inheritance. This can account for both double reduction and preferential pairing.}
#'     \item{\code{"auto"}}{Autopolyploid. Allows only for polysomic inheritance. No double reduction.}
#'     \item{\code{"auto_dr"}}{Autopolyploid, allowing for the effects of double reduction.}
#'     \item{\code{"allo"}}{Allopolyploid. Only complete disomic inheritance is explored.}
#'     \item{\code{"allo_pp"}}{Allopolyploid, allowing for the effects of partial preferential pairing.}
#'     \item{\code{"auto_allo"}}{Only complete disomic and complete polysomic inheritance is studied.}
#'   }
#' @param outlier A logical. Should we allow for outliers (\code{TRUE}) or not (\code{FALSE})?
#' @param ob The default upper bound on the outlier proportion.
#' @param db Should we use the complete equational segregation model (\code{"ces"}) or
#'    the pure random chromatid segregation model (\code{"prcs"}) to determine the upper
#'    bound(s) on the double reduction rate(s). See \code{\link{drbounds}()}
#'    for details.
#' @param ntry The number of times to try the optimization
#'
#' @author David Gerard
#'
#' @examples
#' p1_ploidy <- 4
#' p1 <- 1
#' p2_ploidy <- 8
#' p2 <- 4
#' q <- gf_freq(
#'   p1_g = p1,
#'   p1_ploidy = p1_ploidy,
#'   p1_gamma = 1,
#'   p1_type = "mix",
#'   p2_g = p2,
#'   p2_ploidy = p2_ploidy,
#'   p2_gamma= c(0.2, 0.2, 0.6),
#'   p2_type = "mix",
#'   pi = 0.01)
#' nvec <- c(stats::rmultinom(n = 1, size = 200, prob = q))
#' gl <- simgl(nvec = nvec)
#' seg_lrt(x = nvec, p1_ploidy = p1_ploidy, p2_ploidy = p2_ploidy, p1 = p1, p2 = p2)$p_value
#' seg_lrt(x = gl, p1_ploidy = p1_ploidy, p2_ploidy = p2_ploidy, p1 = p1, p2 = p2)$p_value
#'
#'
#' seg_lrt(x = c(1L, 23L, 51L, 23L, 1L, 1L, 0L), p1_ploidy = 4, p2_ploidy = 8, p1 = 2, p2 = 2)$p_value
#' seg_lrt(x = c(0L, 39L, 115L, 42L, 3L, 1L, 0L), p1_ploidy = 4, p2_ploidy = 8, p1 = 2, p2 = 2)$p_value
#'
#' @export
seg_lrt <- function(
    x,
    p1_ploidy,
    p2_ploidy = p1_ploidy,
    p1 = NULL,
    p2 = NULL,
    model = c("seg", "auto", "auto_dr", "allo", "allo_pp", "auto_allo"),
    outlier = TRUE,
    ob = 0.03,
    db = c("ces", "prcs"),
    ntry = 10) {

  ## Check input ----------
  model <- match.arg(model)
  db <- match.arg(db)

  stopifnot(
    p1_ploidy %% 2 == 0, p1_ploidy > 1,
    p2_ploidy %% 2 == 0, p2_ploidy > 1)
  ## Check if genotype counts or genotype log likelihoods.
  if (is.matrix(x)) {
    stopifnot(ncol(x) == (p1_ploidy + p2_ploidy) / 2 + 1)
    data_type <- "glike"
  } else {
    stopifnot(length(x) == (p1_ploidy + p2_ploidy) / 2 + 1, x >= 0)
    data_type <- "gcount"
  }
  ## Get p1 data type. p1_pos is possible dosages, p1 is genotype log-likelihoods
  if (is.null(p1)) {
    p1 <- rep(0, times = p1_ploidy + 1)
    p1_pos <- 0:p1_ploidy
  } else if (length(p1) == 1) {
    stopifnot(p1 >= 0, p1 <= p1_ploidy, length(p1) == 1)
    p1_pos <- p1
    p1 <- rep(-Inf, times = p1_ploidy + 1)
    p1[p1_pos + 1] <- 0
  } else {
    stopifnot(length(p1) == p1_ploidy + 1)
    p1_pos <- 0:p1_ploidy
  }
  ## Get p2 data type. p2_pos is possible dosages, p2 is genotype log-likelihoods
  if (is.null(p2)) {
    p2 <- rep(0, times = p2_ploidy + 1)
    p2_pos <- 0:p2_ploidy
  } else if (length(p2) == 1) {
    stopifnot(p2 >= 0, p2 <= p2_ploidy, length(p2) == 1)
    p2_pos <- p2
    p2 <- rep(-Inf, times = p2_ploidy + 1)
    p2[p2_pos + 1] <- 0
  } else {
    stopifnot(length(p2) == p2_ploidy + 1)
    p2_pos <- 0:p2_ploidy
  }
  ## check outlier inputs
  stopifnot(length(outlier) == 1, is.logical(outlier))
  stopifnot(length(ob) == 1, ob >= 0, ob <= 1)

  ## Set up gam
  gam <- list()
  gam[[1]] <- list(
    ploidy = p1_ploidy,
    g = NULL,
    alpha = NULL,
    beta = NULL,
    gamma = NULL,
    type = NULL,
    add_dr = NULL)
  gam[[2]] <- list(
    ploidy = p2_ploidy,
    g = NULL,
    alpha = NULL,
    beta = NULL,
    gamma = NULL,
    type = NULL,
    add_dr = NULL)
  gam[[3]] <- list(
    outlier = outlier,
    pi = NULL
  )
  if (model == "seg") {
    gam[[1]]$type <- "mix"
    gam[[1]]$add_dr <- TRUE
    gam[[2]]$type <- "mix"
    gam[[2]]$add_dr <- TRUE
  } else if (model == "auto") {
    gam[[1]]$type <- "polysomic"
    gam[[1]]$add_dr <- FALSE
    gam[[2]]$type <- "polysomic"
    gam[[2]]$add_dr <- FALSE
  } else if (model == "auto_dr") {
    gam[[1]]$type <- "polysomic"
    gam[[1]]$add_dr <- TRUE
    gam[[2]]$type <- "polysomic"
    gam[[2]]$add_dr <- TRUE
  } else if (model == "allo" || model == "allo_pp" || model == "auto_allo") {
    gam[[1]]$type <- "mix"
    gam[[1]]$add_dr <- FALSE
    gam[[2]]$type <- "mix"
    gam[[2]]$add_dr <- FALSE
  }

  ## Alternative optimization ----
  if (data_type == "glike") {
    log_qhat1 <- c(em_li(B = x))
    qhat1 <- exp(log_qhat1)
    l1 <- llike_li(B = x, lpivec = log_qhat1)
  } else if (data_type == "gcount") {
    qhat1 <- x / sum(x)
    l1 <- stats::dmultinom(x = x, prob = x / sum(x), log = TRUE)
  }
  alt_best <- list(
    l1 = l1,
    q1 = qhat1
  )

  ## Null optimization ----
  if (model == "seg" || model == "allo_pp" || model == "auto" || model == "auto_dr") {
    fix_list <- list(
      list(
        alpha = FALSE,
        beta = FALSE,
        gamma = FALSE
      ),
      list(
        alpha = FALSE,
        beta = FALSE,
        gamma = FALSE
      ),
      list(
        pi = FALSE
      )
    )

    if (model == "auto") {
      fix_list[[1]]$alpha <- TRUE
      fix_list[[2]]$alpha <- TRUE
      gam[[1]]$alpha <- rep(x = 0, times = floor(gam[[1]]$ploidy / 4))
      gam[[2]]$alpha <- rep(x = 0, times = floor(gam[[1]]$ploidy / 4))
    }

    null_best <- list()
    null_best$l0_pp <- -Inf
    for (p1_geno in p1_pos) {
      for (p2_geno in p2_pos) {
        gam[[1]]$g <- p1_geno
        gam[[2]]$g <- p2_geno
        nudge <- pkg_env$TOL_small ## bound for starting values
        for (i in seq_len(ntry)) {
          ## Do parent 1
          if (model == "auto_dr") {
            ndr <- floor(gam[[1]]$ploidy / 4)
            gam[[1]]$alpha <- stats::runif(n = ndr, min = rep(x = nudge, times = ndr), max = drbounds(ploidy = gam[[1]]$ploidy, model = db) - nudge)
          } else if (model == "seg") {
            if (p1_geno == 1 || p1_geno == p1_ploidy - 1) {
              gam[[1]]$gamma <- 1
              gam[[1]]$beta <- stats::runif(n = 1, min = nudge, max = beta_bounds(ploidy = p1_ploidy) - nudge)
            } else {
              nmix <- n_pp_mix(g = p1_geno, ploidy = p1_ploidy)
              gam[[1]]$gamma <- stats::runif(n = nmix)
              gam[[1]]$gamma <- gam[[1]]$gamma / sum(gam[[1]]$gamma)
            }
          } else if (model == "allo_pp") {
            if (p1_geno == 1 || p1_geno == p1_ploidy - 1) {
              gam[[1]]$gamma <- 1
            } else {
              nmix <- n_pp_mix(g = p1_geno, ploidy = p1_ploidy)
              gam[[1]]$gamma <- stats::runif(n = nmix)
              gam[[1]]$gamma <- gam[[1]]$gamma / sum(gam[[1]]$gamma)
            }
          }

          ## do parent 2
          if (model == "auto_dr") {
            ndr <- floor(gam[[2]]$ploidy / 4)
            gam[[2]]$alpha <- stats::runif(n = ndr, min = rep(x = nudge, times = ndr), max = drbounds(ploidy = gam[[2]]$ploidy, model = db) - nudge)
          } else if (model == "seg") {
            if (p2_geno == 1 || p2_geno == p2_ploidy - 1) {
              gam[[2]]$gamma <- 1
              gam[[2]]$beta <- stats::runif(n = 1, min = nudge, max = beta_bounds(ploidy = p2_ploidy) - nudge)
            } else {
              nmix <- n_pp_mix(g = p2_geno, ploidy = p2_ploidy)
              gam[[2]]$gamma <- stats::runif(n = nmix)
              gam[[2]]$gamma <- gam[[2]]$gamma / sum(gam[[2]]$gamma)
            }
          } else if (model == "allo_pp") {
            if (p2_geno == 1 || p2_geno == p2_ploidy - 1) {
              gam[[2]]$gamma <- 1
            } else {
              nmix <- n_pp_mix(g = p2_geno, ploidy = p2_ploidy)
              gam[[2]]$gamma <- stats::runif(n = nmix)
              gam[[2]]$gamma <- gam[[2]]$gamma / sum(gam[[2]]$gamma)
            }
          }

          if (gam[[3]]$outlier) {
            gam[[3]]$pi <- stats::runif(n = 1, min = 0, max = ob)
          }

          ret <- gam_to_par(gam = gam, fix_list = fix_list, db = db, ob = ob)
          oout <- stats::optim(
            par = ret$par,
            fn = ifelse(data_type == "glike", ll_gl, ll_g),
            method = ifelse(length(ret$par) == 1, "Brent", "L-BFGS-B"),
            rule = ret$rule,
            x = x,
            control = list(fnscale = -1),
            lower = ret$lower,
            upper = ret$upper)

          if (oout$value + p1[[p1_geno + 1]] + p2[[p2_geno + 1]] > null_best$l0_pp) {
            null_best$l0_pp <- oout$value + p1[[p1_geno + 1]] + p2[[p2_geno + 1]]
            null_best$l0 <- oout$value
            null_best$q0 <- par_to_gf(par = oout$par, rule = ret$rule) ## ML genotype frequency
            null_best$gam <- par_to_gam(par = oout$par, rule = ret$rule) ## MLE's
            ## find df0
            if (model == "seg" || model == "allo_pp") {
              null_best$df0 <- gam_to_df_seg(gam = null_best$gam, db = db, ob = ob)
            } else if (model == "auto") {
              if (outlier) {
                if (null_best$gam[[3]]$pi > pkg_env$TOL_big && null_best$gam[[3]]$pi < ob - pkg_env$TOL_big) {
                  null_best$df0 <- 1
                } else {
                  null_best$df0 <- 0
                }
              } else {
                null_best$df0 <- 0
              }
            } else if (model == "auto_dr") {
                null_best$df0 <- gam_to_df_auto_dr(gam = null_best$gam, db = db, ob = ob)
            }
          }
        }
      }
    }
  } else if (model == "allo" || model == "auto_allo") {

  } else {
    stop("seg_lrt: how did you get here?")
  }

  ## subtract off if both zero under null and alt, then subtract off one.
  # alt_best$df1 <- length(alt_best$q1) - 1
  alt_best$df1 <- length(alt_best$q1) - sum((alt_best$q1 < pkg_env$TOL_big) & (null_best$q0 < pkg_env$TOL_big)) - 1

  ## Run test
  ret <- list(
    null = null_best,
    alt = alt_best,
    stat = 2 * (alt_best$l1 - null_best$l0),
    df = max(alt_best$df1 - null_best$df0, 1)
  )
  ret$p_value <- stats::pchisq(q = ret$stat, df = ret$df, lower.tail = FALSE)

  return(ret)
}


#' Calculate number of parameters under the null.
#'
#' Does so by counting number on boundary and then passing those to
#' gam_to_nparam_seg to get identified dimesion.
#'
#' @inheritParams gam_to_par
#'
#' This only works for model = "seg" and model = "allo_pp"
#'
#' @author David Gerard
#'
#' @noRd
gam_to_df_seg <- function(gam, fix_list = NULL, db = c("ces", "prcs"), ob = 0.03) {
  stopifnot(is.null(gam[[1]]$alpha) || is.null(gam[[1]]$beta))
  stopifnot(is.null(gam[[1]]$alpha) || is.null(gam[[1]]$gamma))
  stopifnot(is.null(gam[[2]]$alpha) || is.null(gam[[2]]$beta))
  stopifnot(is.null(gam[[2]]$alpha) || is.null(gam[[2]]$gamma))
  TOL <- pkg_env$TOL_big
  db <- match.arg(db)
  if (is.null(fix_list)) {
    fix_list <- list(
      list(
        alpha = FALSE,
        beta = FALSE,
        gamma = FALSE
      ),
      list(
        alpha = FALSE,
        beta = FALSE,
        gamma = FALSE
      ),
      list(
        pi = FALSE
      )
    )
  }

  df <- 0

  if (gam[[1]]$g == 0 || gam[[1]]$g == gam[[1]]$ploidy) {
    ## do nothing
  } else if (gam[[1]]$g == 1 || gam[[1]]$g == gam[[1]]$ploidy - 1) {
    if (!is.null(gam[[1]]$beta) && !fix_list[[1]]$beta) {
      if (gam[[1]]$beta > TOL && gam[[1]]$beta < beta_bounds(ploidy = gam[[1]]$ploidy, model = db) - TOL) {
        df <- df + 1
      }
    } else if (!is.null(gam[[1]]$alpha) && !fix_list[[1]]$alpha) {
      beta <- sum(seq_along(gam[[1]]$alpha) * gam[[1]]$alpha) / gam[[1]]$ploidy
      if (beta > TOL && beta < beta_bounds(ploidy = gam[[1]]$ploidy, model = db) - TOL) {
        df <- df + 1
      }
    }
  } else {
    if (!is.null(gam[[1]]$gamma) && !fix_list[[1]]$gamma) {
      df <- df + sum(gam[[1]]$gamma > TOL) - 1
    } else if (!is.null(gam[[1]]$alpha && !fix_list[[1]]$alpha)) { ## only part not 100% sure about for very large ploidies. Double check.
      df <- df + sum((gam[[1]]$alpha > TOL) & (gam[[1]]$alpha < drbounds(ploidy = gam[[1]]$ploidy, model = db)))
    }
  }

  if (gam[[2]]$g == 0 || gam[[2]]$g == gam[[2]]$ploidy) {
    ## do nothing
  } else if (gam[[2]]$g == 1 || gam[[2]]$g == gam[[2]]$ploidy - 1) {
    if (!is.null(gam[[2]]$beta) && !fix_list[[1]]$beta) {
      if (gam[[2]]$beta > TOL && gam[[2]]$beta < beta_bounds(ploidy = gam[[2]]$ploidy, model = db) - TOL) {
        df <- df + 1
      }
    } else if (!is.null(gam[[2]]$alpha) && !fix_list[[1]]$alpha) {
      beta <- sum(seq_along(gam[[2]]$alpha) * gam[[2]]$alpha) / gam[[2]]$ploidy
      if (beta > TOL && beta < beta_bounds(ploidy = gam[[2]]$ploidy, model = db) - TOL) {
        df <- df + 1
      }
    }
  } else {
    if (!is.null(gam[[2]]$gamma) && !fix_list[[1]]$gamma) {
      df <- df + sum(gam[[2]]$gamma > TOL) - 1
    } else if (!is.null(gam[[2]]$alpha && !fix_list[[1]]$alpha)) {
      df <- df + sum((gam[[2]]$alpha > TOL) & (gam[[2]]$alpha < drbounds(ploidy = gam[[2]]$ploidy, model = db)))
    }
  }

  ## Adjust for special case of two gammas, both in interior, neither fixed
  ## or two betas, both interior, neither fixed.
  ## overwrite df in this case
  if (!fix_list[[1]]$gamma && !fix_list[[2]]$gamma && !is.null(gam[[1]]$gamma) && !is.null(gam[[2]]$gamma) && gam[[1]]$g > 1 && gam[[1]]$g < gam[[1]]$ploidy - 1 && gam[[2]]$g > 1 && gam[[2]]$g < gam[[2]]$ploidy - 1) {
    if (all(gam[[1]]$gamma > TOL) && all(gam[[2]]$gamma > TOL)) {
      df <- gam_to_nparam_seg(gam = gam)
    }
  }
  if (!fix_list[[1]]$beta && !fix_list[[2]]$beta && !is.null(gam[[1]]$beta) && !is.null(gam[[2]]$beta) && (gam[[1]]$g == 1 || gam[[1]]$g == gam[[1]]$ploidy - 1) && (gam[[2]]$g == 1 || gam[[2]]$g == gam[[2]]$ploidy - 1)) {
    if (gam[[1]]$beta > TOL && gam[[1]]$beta < beta_bounds(ploidy = gam[[1]]$ploidy, model = db) - TOL && gam[[2]]$beta > TOL && gam[[2]]$beta < beta_bounds(ploidy = gam[[2]]$ploidy, model = db) - TOL) {
      df <- gam_to_nparam_seg(gam = gam)
    }
  }


  if (gam[[3]]$outlier && !fix_list[[3]]$pi && !is.null(gam[[3]]$pi)) {
    if (gam[[3]]$pi > TOL && gam[[3]]$pi < ob - TOL) {
      df <- df + 1
    }
  }

  return(df)
}

#' Calculates the true dimension of the parameter space
#'
#' Each parent should have non-null gamma parameters. So no simplex
#' and nullplex parents.
#'
#' This only works for model = "seg" or model = "allo_pp"
#'
#' The parameterspace might not be identified. This calculates the jacobian of the
#' genotype frequencies (dependent variable) as a function of the null parameters
#' (independent variable). The rank of the jacobian is the true number of parameters.
#'
#' @inheritParams gam_to_par
#'
#' @author David Gerard
#'
#' @noRd
gam_to_nparam_seg <- function(gam) {
  if (gam[[1]]$g > 1 && gam[[1]]$g < gam[[1]]$ploidy - 1 && gam[[2]]$g > 1 && gam[[2]]$g < gam[[2]]$ploidy - 1) {
    fix_list <- list(
      list(
        alpha = TRUE,
        beta = TRUE,
        gamma = FALSE
      ),
      list(
        alpha = TRUE,
        beta = TRUE,
        gamma = FALSE
      ),
      list(
        pi = TRUE
      )
    )
  } else if ((gam[[1]]$g == 1 || gam[[1]]$g == gam[[1]]$ploidy - 1) && (gam[[2]]$g == 1 || gam[[2]]$g == gam[[2]]$ploidy - 1)) {
    fix_list <- list(
      list(
        alpha = TRUE,
        beta = FALSE,
        gamma = TRUE
      ),
      list(
        alpha = TRUE,
        beta = FALSE,
        gamma = TRUE
      ),
      list(
        pi = TRUE
      )
    )
  } else {
    stop("gam_to_nparam_seg: only both simplex or both multiplex allowed?")
  }

  ret <- gam_to_par(gam = gam, fix_list = fix_list)

  env <- new.env()
  env$par <- ret$par
  env$rule <- ret$rule
  dout <- stats::numericDeriv(expr = quote(par_to_gf(par = par, rule = rule)), theta = "par", rho = env)
  return(sum(svd(attr(dout, "gradient"))$d > pkg_env$TOL_big))
}

#' Number of null parameters under auto_dr model
#'
#' Counts the dimenson of the inner parameters
#'
#' @inheritParams gam_to_par
#'
#' @author David Gerard
#'
#' @noRd
gam_to_df_auto_dr <- function(gam,  db = c("ces", "prcs"), ob = 0.03) {
  db <- match.arg(db)
  TOL <- pkg_env$TOL_big

  ret <- gam_to_par(gam = gam, db = db, ob = ob)
  onbound <- (ret$par < ret$lower + TOL) | (ret$par > ret$upper - TOL)
  if (sum(!onbound) == 0) {
    return(0)
  }
  par_inner <- ret$par[!onbound]
  par_edge <-ret$par[onbound]
  fn <- function(par_inner, par_edge, onbound, rule) {
    par <- rep(NA_real_, times = length(onbound))
    par[!onbound] <- par_inner
    par[onbound] <- par_edge
    par_to_gf(par = par, rule = rule)
  }

  env <- new.env()
  env$par_inner <- par_inner
  env$par_edge <- par_edge
  env$onbound <- onbound
  env$rule <- ret$rule
  dout <- stats::numericDeriv(expr = quote(fn(par_inner = par_inner, par_edge = par_edge, onbound = onbound, rule = rule)), theta = "par_inner", rho = env)
  nparam <- sum(svd(attr(dout, "gradient"))$d > TOL)
  return(nparam)
}

