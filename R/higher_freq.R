logit <- function (x) {
    log(x) - log(1 - x)
}

expit <- function (x) {
    1/(1 + exp(-x))
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


f1_gf_pp <- function(gamma1, gamma2, g1, g2) {

}
