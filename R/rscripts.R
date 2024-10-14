#' Tetraploid gamete frequencies of gametes when one parent's genotype is known
#'
#' This is under the two parameter model.
#'
#' @param alpha The double reduction rate
#' @param xi The preferential pairing parameter
#' @param ell The parental genotype
#'
#' @return The gamete genotype frequencies
#'
#' @author Mira Thakkar and David Gerard
#'
#' @examples
#' alpha <- 1/6
#' xi <- 1/3
#' pvec_tet_2(alpha = alpha, xi = xi, ell = 0)
#' pvec_tet_2(alpha = alpha, xi = xi, ell = 1)
#' pvec_tet_2(alpha = alpha, xi = xi, ell = 2)
#' pvec_tet_2(alpha = alpha, xi = xi, ell = 3)
#' pvec_tet_2(alpha = alpha, xi = xi, ell = 4)
#'
#' @export
pvec_tet_2 <- function(alpha, xi, ell) {
  if (ell > 4 | ell < 0 | is.na(ell)){
    stop("Invalid input")
  } else if(ell == 0){
    pv <- c(
      1,
      0,
      0
    )
  } else if(ell == 1){
    pv <- c(
      0.5 + 0.25 * alpha,
      0.5 - 0.5 * alpha,
      0.25 * alpha
    )
  } else if(ell == 2){
    pv <- c(
      0.5 * alpha + 0.25 * (1 - alpha) * (1 - xi),
      0.5 * (1 - alpha) * (1 + xi),
      0.5 * alpha + 0.25 * (1 - alpha) * (1 - xi)
    )
  } else if(ell == 3){
    pv <- c(
      0.25 * alpha,
      0.5 - 0.5 * alpha,
      0.5 + 0.25 * alpha
    )
  } else if(ell == 4){
    pv <- c(
      0,
      0,
      1
    )
  }
  return(pv)
}

#' Tetraploid gamete frequencies of gametes when one parent's genotype is known
#'
#' This is under the three parameter model.
#'
#' @param tau Probability of quadrivalent formation
#' @param beta Probability of double reduction given quadrivalent formation
#' @param gamma Probability of AA/aa pairing given bivalent formation
#' @param ell The parent genotype
#'
#' @return The gamete genotype frequencies
#'
#' @author David Gerard
#'
#' @examples
#' pvec_tet_3(tau = 0.5, beta = 0.1, gamma = 0.5, ell = 2)
#'
#' @export
pvec_tet_3 <- function(tau, beta, gamma, ell) {
  parvec <- three_to_two(tau = tau, beta = beta, gamma = gamma)
  gf <- pvec_tet_2(alpha = parvec[["alpha"]], xi = parvec[["xi"]], ell = ell)
  return(gf)
}


#' Calculates offspring genotype frequencies under the two-parameter model.
#'
#' @param alpha The double reduction rate
#' @param xi1 The preferential pairing parameter of the first parent.
#' @param xi2 The preferential pairing parameter of the second parent.
#' @param p1 The first parent's genotype
#' @param p2 The second parent's genotype
#'
#' @return Offspring genotype frequencies
#'
#' @author Mira Thakkar
#'
#' @examples
#' alpha <- 1/6
#' xi1 <- 1/3
#' xi2 <- 1/3
#' p1 <- 2
#' p2 <- 3
#' offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = p1, p2 = p2)
#'
#' @export
offspring_gf_2 <- function(alpha, xi1, xi2 = xi1, p1, p2){

  pvec1 <- pvec_tet_2(alpha = alpha, xi = xi1, ell = p1)
  pvec2 <- pvec_tet_2(alpha = alpha, xi = xi2, ell = p2)

  qvec <- stats::convolve(pvec1, rev(pvec2), type = "open")

  stopifnot(qvec > -1e-06)
  stopifnot(abs(sum(qvec) - 1) < 1e-06)
  TOL <- sqrt(.Machine$double.eps)
  qvec[qvec < TOL] <- 0
  qvec <- qvec / sum(qvec)

  return(qvec)
}

#' Calculates offspring genotype frequencies under the three-parameter model.
#'
#' @inheritParams pvec_tet_3
#' @inheritParams offspring_gf_2
#' @param gamma1 Probability of AA_aa pairing in parent 1
#' @param gamma2 Probability of AA_aa pairing in parent 2
#'
#' @return Offspring genotype frequencies
#'
#' @author David Gerard
#'
#' @examples
#' offspring_gf_3(
#'   tau = 1/2,
#'   beta = 1/6,
#'   gamma1 = 1/3,
#'   gamma2 = 1/3,
#'   p1 = 1,
#'   p2 = 2)
#'
#' @export
offspring_gf_3 <- function(tau, beta, gamma1, gamma2 = gamma1, p1, p2) {
  parvec1 <- three_to_two(tau = tau, beta = beta, gamma = gamma1)
  parvec2 <- three_to_two(tau = tau, beta = beta, gamma = gamma2)
  stopifnot(parvec1[["alpha"]] == parvec2[["alpha"]])
  gf <- offspring_gf_2(
    alpha = parvec1[["alpha"]],
    xi1 = parvec1[["xi"]],
    xi2 = parvec2[["xi"]],
    p1 = p1,
    p2 = p2)
  return(gf)
}

#' Convert from three parameters to two parameters
#'
#' @inheritParams pvec_tet_3
#'
#' @return A vector of length two. The first is the double reduction rate
#'    (\code{alpha}), and the second is the preferential pairing parameter
#'    (\code{xi}).
#'
#' @author David Gerard
#'
#' @examples
#' three_to_two(tau = 0.1, beta = 1/6, gamma = 1/4)
#'
#' @export
three_to_two <- function(tau, beta, gamma) {
  alpha <- tau * beta
  eta <- (1 - beta) * tau / ((1 - beta) * tau + (1 - tau))
  xi <- eta / 3 + (1 - eta) * gamma
  return(c(alpha = alpha, xi = xi))
}


#' Simulates genotypes given genotype frequencies.
#'
#' Takes as input the offspring
#' genotype frequencies and a sample size and returns simulated genotypes.
#'
#' @param gf Vector of offspring genotype frequencies
#' @param n Sample size
#'
#' @return Simulated genotypes
#'
#' @author Mira Thakkar
#'
#' @examples
#' gf <- offspring_gf_2(alpha = 1/6, xi1 = 1/3, xi2 = 1/3, p1 = 2, p2 = 3)
#' offspring_geno(gf = gf, n = 10)
#'
#' @export
offspring_geno <- function(gf, n){
  sim_gen <- c(stats::rmultinom(n = 1, size = n, prob = gf))
  return(sim_gen)
}

#' Converts genotype counts to genotype vectors.
#'
#' @param gcount The vector of genotype counts.
#'
#' @author David Gerard
#'
#' @export
#'
#' @seealso [gvec_to_gcount()]
#'
#' @examples
#' gcount <- c(1, 2, 3, 0, 5)
#' gcount_to_gvec(gcount = gcount)
gcount_to_gvec <- function(gcount) {
  unlist(mapply(FUN = rep, x = seq_along(gcount) - 1, each = gcount))
}

#' Inverse function of \code{\link{gcount_to_gvec}()}.
#'
#' @param gvec The vector of genotypes. \code{gvec[i]} is the genotype
#'     for individual i.
#' @param ploidy The ploidy of the species.
#'
#' @author David Gerard
#'
#' @export
#'
#' @seealso [gcount_to_gvec()]
#'
#' @examples
#' gvec <- c(1, 2, 3, 2, 3, 1, 4, 0, 1, 0, 0, 1, 0, 0)
#' gvec_to_gcount(gvec = gvec)
#'
gvec_to_gcount <- function(gvec, ploidy = 4) {
  x <- table(factor(gvec, levels = 0:ploidy))
  attributes(x) <- NULL
  return(x)
}

#' Generate genotype likelihoods from offspring genotypes.
#'
#' Takes as input (i) the parent genotypes,
#' (ii) the offspring genotype freq, (iii) sequencing error rate, (iv) read
#' depth, (v) bias, (vi) overdispersion and returns genotype likelihoods.
#'
#' @param genovec Offspring genotypes. \code{genovec[i]} is the dosage for individual i.
#' @param p1_geno Parent 1 genotype
#' @param p2_geno Parent 2 genotype
#' @param ploidy Ploidy
#' @param rd Read depth. Lower is more uncertain.
#' @param seq Sequencing error rate. Higher means more uncertain.
#' @param bias Bias. 1 means no bias.
#' @param od Overdispersion. Typical value is like 0.01. Higher means more uncertain.
#'
#' @return Genotype likelihoods
#'
#' @author Mira Thakkar
#'
#' @export
po_gl <- function(genovec, ploidy, p1_geno = NULL, p2_geno = NULL, rd = 10, seq = 0.01, bias = 1, od = 0.01) {
  n <- length(genovec)
  stopifnot(genovec <= ploidy)
  sizevec <- rep(rd, length.out = n)
  refvec <- updog::rflexdog(sizevec = sizevec, geno = genovec, ploidy = ploidy, seq = seq, bias = bias, od = od)
  if (!is.null(p1_geno)) {
    p1ref <- updog::rflexdog(sizevec = rd, geno = p1_geno, ploidy = ploidy, seq = seq, bias = bias, od = od)
    p1rd <- rd
  } else {
    p1ref <- NULL
    p1rd <- NULL
  }

  if (!is.null(p2_geno)) {
    p2ref <- updog::rflexdog(sizevec = rd, geno = p2_geno, ploidy = ploidy, seq = seq, bias = bias, od = od)
    p2rd <- rd
  } else {
    p2ref <- NULL
    p2rd <- NULL
  }

  fout <- updog::flexdog_full(refvec = refvec,
                              sizevec = sizevec,
                              ploidy = ploidy,
                              model = "f1pp",
                              seq = seq,
                              bias = bias,
                              od = od,
                              update_bias = FALSE,
                              update_seq = FALSE,
                              update_od = FALSE,
                              p1ref = p1ref,
                              p1size = p1rd,
                              p2ref = p2ref,
                              p2size = p2rd)

  return(fout)
}

#' Simulate genotype likelihoods of F1 individuals.
#'
#' @param n Sample size.
#' @param g1 The first parent's genotype.
#' @param g2 The second parent's genotype.
#' @param rd The read depth.
#' @param alpha The double reduction rate.
#' @param xi1 The first parent's preferential pairing parameter.
#' @param xi2 The second parent's preferential pairing parameter.
#'
#' @return The matrix of offspring genotype log-likelihoods.
#'
#' @examples
#' set.seed(1)
#' simf1gl(n = 10, g1 = 1, g2 = 2)
#'
#' @author David Gerard
#'
#' @export
simf1gl <- function(n, g1, g2, rd = 10, alpha = 0, xi1 = 1/3, xi2 = 1/3) {
  gcount <- simf1g(n = n, g1 = g1, g2 = g2, alpha = alpha, xi1 = xi1, xi2 = xi2)
  gvec <- gcount_to_gvec(gcount = gcount)
  fout <- po_gl(genovec = gvec, p1_geno = g1, p2_geno = g2, ploidy = 4, rd = rd)
  rownames(fout$genologlike) <- paste0("F", seq_len(nrow(fout$genologlike)))
  colnames(fout$genologlike) <- 0:4
  return(fout$genologlike)
}

#' Simulate genotype counts from F1 individuals
#'
#' @inheritParams simf1gl
#'
#' @return A vector of counts, where element \code{i} is the number of
#'     simulated individuals with genotype \code{i-1}.
#'
#' @author David Gerard
#'
#' @examples
#' set.seed(1)
#' simf1g(n = 10, g1 = 1, g2 = 2)
#'
#' @export
simf1g <- function(n, g1, g2, alpha = 0, xi1 = 1/3, xi2 = 1/3) {
  gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = g1, p2 = g2)
  gcount <- offspring_geno(gf = gf, n = n)
  names(gcount) <- 0:4
  return(gcount)
}

#' Simulate genotype likelihoods from genotype counts
#'
#' Provide a vector of genotype counts and this will return a matrix of
#' genotype log-likelihoods.
#'
#' @param nvec A vector of counts. \code{nvec[k]} is the number of folks with a genotype of k-1.
#' @inheritParams po_gl
#'
#' @return A matrix of genotype log-likelihoods. The rows index the
#'      individuals and the columns index the genotypes. This is natural
#'      log (base e).
#'
#' @author David Gerard
#'
#' @export
simgl <- function(nvec, rd = 10, seq = 0.01, bias = 1, od = 0.01) {
  ploidy <- length(nvec) - 1
  genovec <- gcount_to_gvec(gcount = nvec)
  ret <- po_gl(genovec = genovec, ploidy = ploidy, p1_geno = NULL, p2_geno = NULL, seq = seq, rd = rd, bias = bias, od = od)
  return(ret$genologlike)
}
