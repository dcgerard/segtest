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
    mat_list[[i + 1]] <- cbind(i, all_multinom(n - i, k -
                                                 1))
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
    if (all(abs(alpha) < sqrt(.Machine$double.eps))) {
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
  TOL <- sqrt(.Machine$double.eps)
  stopifnot(
    length(gamma) == n_pp_mix(g = g, ploidy = ploidy),
    abs(sum(gamma) - 1) < TOL,
    gamma >= 0)
  plist <- seg[seg$ploidy == ploidy & (seg$mode == "disomic" | seg$mode == "both") & seg$g == g, ]$p
  pvec <- rep(0, length.out = ploidy / 2 + 1)
  for (i in seq_along(plist)) {
    pvec <- pvec + plist[[i]] * gamma[[i]]
  }
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
  TOL <- sqrt(.Machine$double.eps)
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

#' Convert parameterization from something optim() can use to something
#' gamfreq() can use.
#'
#' @param rule A list of length three. The first element gives the model
#'     of parent 1. The second element gives the model of parent 2.
#'     The third element is a logiical on whether we add an outlier or now.
#'     \describe {
#'       \item{ploidy}{The parent's ploidy.}
#'       \item{g}{The parent dosage.}
#'       \item{type}{Either "mix", "mix_dr", or "polysomic".}
#'       \item{outlier}{(only in third element) A logical on whether outliers are modeled.}
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
    alpha = NULL,
    beta = NULL,
    gamma = NULL,
    type = NULL,
    add_dr = NULL)
  gam[[2]] <- list(
    ploidy = rule[[2]]$ploidy,
    g = rule[[2]]$g,
    alpha = NULL,
    beta = NULL,
    gamma = NULL,
    type = NULL,
    add_dr = NULL)
  gam[[3]] <- list(
    outlier = rule[[3]]$outlier,
    pi = NULL
  )

  ## Do parent 1 ----
  cindex <- 1 ## keeps track of where we are in par
  if (rule[[1]]$type == "polysomic") {
    gam[[1]]$type <- "polysomic"
    gam[[1]]$add_dr <- FALSE
    if (rule[[1]]$g != 0 && rule[[1]]$g != rule[[1]]$ploidy) {
      ndr <- floor(rule[[1]]$ploidy / 4)
      gam[[1]]$alpha <- par[cindex:(cindex + ndr - 1)]
      cindex <- cindex + ndr
    }
  } else if (rule[[1]]$type == "mix_dr" && (rule[[1]]$g == 1 || rule[[1]]$g == rule[[1]]$ploidy - 1)) {
    gam[[1]]$type <- "mix"
    gam[[1]]$add_dr <- TRUE
    gam[[1]]$beta <- par[[cindex]]
    cindex <- cindex + 1
  } else if (rule[[1]]$type == "mix" && (rule[[1]]$g == 1 || rule[[1]]$g == rule[[1]]$ploidy - 1)) {
    gam[[1]]$type <- "mix"
    gam[[1]]$gamma <- 1
    gam[[1]]$add_dr <- FALSE
  } else if ((rule[[1]]$type == "mix" || rule[[1]]$type == "mix_dr")) {
    gam[[1]]$type <- "mix"
    gam[[1]]$add_dr <- FALSE
    if (rule[[1]]$g != 0 && rule[[1]]$g != rule[[1]]$ploidy) {
      npp <- n_pp_mix(g = rule[[1]]$g, ploidy = rule[[1]]$ploidy)
      gam[[1]]$gamma <- real_to_simplex(par[cindex:(cindex + npp - 2)])
      cindex <- cindex + npp - 1
    }
  } else {
    stop("par_to_gam 1: how did you get here?")
  }

  ## Do parent 2 ----
  if (rule[[2]]$type == "polysomic") {
    gam[[2]]$type <- "polysomic"
    gam[[2]]$add_dr <- FALSE
    if (rule[[2]]$g != 0 && rule[[2]]$g != rule[[2]]$ploidy) {
      ndr <- floor(rule[[2]]$ploidy / 4)
      gam[[2]]$alpha <- par[cindex:(cindex + ndr - 1)]
      cindex <- cindex + ndr
    }
  } else if (rule[[2]]$type == "mix_dr" && (rule[[2]]$g == 1 || rule[[2]]$g == rule[[2]]$ploidy - 1)) {
    gam[[2]]$type <- "mix"
    gam[[2]]$add_dr <- TRUE
    gam[[2]]$beta <- par[[cindex]]
    cindex <- cindex + 1
  } else if (rule[[2]]$type == "mix" && (rule[[2]]$g == 1 || rule[[2]]$g == rule[[2]]$ploidy - 1)) {
    gam[[2]]$type <- "mix"
    gam[[2]]$gamma <- 1
    gam[[2]]$add_dr <- FALSE
  } else if ((rule[[2]]$type == "mix" || rule[[2]]$type == "mix_dr")) {
    gam[[2]]$type <- "mix"
    gam[[2]]$add_dr <- FALSE
    if (rule[[2]]$g != 0 && rule[[2]]$g != rule[[2]]$ploidy) {
      npp <- n_pp_mix(g = rule[[2]]$g, ploidy = rule[[2]]$ploidy)
      gam[[2]]$gamma <- real_to_simplex(par[cindex:(cindex + npp - 2)])
      cindex <- cindex + npp - 1
    }
  } else {
    stop("par_to_gam 2: how did you get here?")
  }

  if (rule[[3]]$outlier) {
    gam[[3]]$pi <- par[[cindex]]
    cindex <- cindex + 1
  }

  return(gam)
}

#' Par to genotype frequencies
#'
#' @inheritParams par_to_gam
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
par_to_gf <- function(par, rule) {
  TOL <- sqrt(.Machine$double.eps)
  gampar <- par_to_gam(par = par, rule = rule)
  p1 <- do.call(what = gamfreq, args = gampar[[1]])
  p2 <- do.call(what = gamfreq, args = gampar[[2]])
  q <- stats::convolve(x = p1, y = rev(p2), type = "open")
  q[q < TOL] <- 0 ## for -1e-16
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
#'     gamfreq of parent 2. The third contains the outlier proportion.
#'     These parameters in the parent lists are \code{ploidy}, \code{g}, \code{gamma},
#'     \code{beta}, \code{alpha}, \code{type}, and \code{add_dr}.
#'     See \code{\link{gamfreq}()} for details on these parameters. The third
#'     is a list with elements \code{outlier} (a logical on whether
#'     outliers are modeled) and \code{pi} (the outlier proportion).
#'
#' @return A list of length 2, containing \code{par} and \code{rule}. See
#'     \code{\link{par_to_gam}()} for a description of these elements.
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
#'   add_dr = NULL)
#' gam[[2]] <- list(
#'   ploidy = 8,
#'   g = 4,
#'   alpha = NULL,
#'   beta = NULL,
#'   gamma = c(0.1, 0.9),
#'   type = "mix",
#'   add_dr = NULL)
#' gam[[3]] <- list(
#'   outlier = TRUE,
#'   pi = 0.03
#' )
gam_to_par <- function(gam) {
  if (gam[[1]]$type == "polysomic") {

  }

}





