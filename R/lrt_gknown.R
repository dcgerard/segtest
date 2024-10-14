#####################
## Likelihood methods when genotypes are known
#####################

#' Likelihood ratio test for segregation distortion with known genotypes
#'
#' This will run a likelihood ratio test using the genotypes of an F1 population
#' of tetraploids for the null of Mendelian segregation (accounting for double
#' reduction and preferential pairing) against the alternative of
#' segregation distortion. This is when the genotypes are assumed known.
#'
#' @section Impossible genotypes:
#' Some offspring genotype combinations are impossible given the parental
#' genotypes. If these impossible genotypes combinations show up, we return a
#' p-value of 0, a log-likelihood ratio statistic of Infinity, and missing
#' values for all other return items. The impossible genotypes are:
#' \describe{
#'   \item{\code{g1 = 0 && g2 = 0}}{Only offspring genotypes of 0 are possible.}
#'   \item{\code{g1 = 4 && g2 = 4}}{Only offspring genotypes of 4 are possible.}
#'   \item{\code{g1 = 0 && g2 = 4 || g1 == 4 && g2 == 0}}{Only offspring genotypes of 2 are possible.}
#'   \item{\code{g1 = 0 && g2 \%in\% c(1, 2, 3) || g1 = \%in\% c(1, 2, 3) && g2 == 0}}{Only offspring genotypes of 0, 1, and 2 are possible.}
#'   \item{\code{g1 = 4 && g2 \%in\% c(1, 2, 3) || g1 = \%in\% c(1, 2, 3) && g2 == 4}}{Only offspring genotypes of 2, 3, and 4 are possible.}
#' }
#'
#' @section Unidentified parameters:
#' When \code{g1 = 2} or \code{g2 = 2} (or both), the model is not identified
#' and those estimates (\code{alpha}, \code{xi1}, and \code{xi2}) are
#' meaningless. Do NOT interpret them.
#'
#' The estimate of \code{alpha} (double reduction rate) IS identified as
#' long as at least one parent is simplex, and no parent is duplex.
#' However, the estimates of the double reduction rate have extremely high
#' variance.
#'
#' @param x A vector of genotype counts. \code{x[i]} is the number of
#'    offspring with genotype \code{i-1}.
#' @param g1 The genotype of parent 1.
#' @param g2 The genotype of parent 2.
#' @param drbound The maximum rate of double reduction. A default of 1/6
#'     is provided, which is the rate under the complete equational
#'     segregation model of meiosis.
#' @param pp A logical. Should we account for preferential pairing
#'     (\code{TRUE}) or not (\code{FALSE})?
#' @param dr A logical. Should we account for double reduction
#'     (\code{TRUE}) or not (\code{FALSE})?
#' @param alpha If \code{dr = FALSE}, this is the known rate of double
#'     reduction.
#' @param xi1 If \code{pp = FALSE}, this is the known preferential pairing
#'     parameter of parent 1.
#' @param xi2 If \code{pp = FALSE}, this is the known preferential pairing
#'     parameter of parent 2.
#'
#' @return A list with the following elements
#' \describe{
#'   \item{\code{statistic}}{The log-likelihood ratio test statistic.}
#'   \item{\code{df}}{The degrees of freedom.}
#'   \item{\code{p_value}}{The p-value.}
#'   \item{\code{alpha}}{The estimated double reduction rate.}
#'   \item{\code{xi1}}{The estimated preferential pairing parameter of parent 1.}
#'   \item{\code{xi2}}{The estimated preferential pairing parameter of parent 2.}
#' }
#'
#' @author David Gerard
#'
#' @examples
#' set.seed(100)
#' gf <- offspring_gf_2(alpha = 1/12, xi1 = 0.2, xi2 = 0.6, p1 = 1, p2 = 0)
#' x <- offspring_geno(gf = gf, n = 100)
#' lrt_men_g4(x = x, g1 = 1, g2 = 0)
#'
#' @export
lrt_men_g4 <- function(
    x,
    g1,
    g2,
    drbound = 1/6,
    pp = TRUE,
    dr = TRUE,
    alpha = 0,
    xi1 = 1/3,
    xi2 = 1/3) {
  ## check input
  stopifnot(length(x) == 5,
            x >= 0,
            g1 >= 0,
            g1 <= 4,
            g2 >= 0,
            g2 <= 4,
            drbound >= 0,
            drbound <= 1,
            alpha >= 0,
            alpha <= 1,
            xi1 >= 0,
            xi1 <= 1,
            xi2 >= 0,
            xi2 <= 1,
            is.logical(pp),
            is.logical(dr),
            length(g1) == 1,
            length(g2) == 1,
            length(drbound) == 1,
            length(pp) == 1,
            length(dr) == 1,
            length(alpha) == 1,
            length(xi1) == 1,
            length(xi2) == 1)

  TOL <- sqrt(.Machine$double.eps)
  if (is_impossible(x = x, g1 = g1, g2 = g2, dr = dr || alpha > TOL)) {
    ## Check if data are possible under null
    ret <- list(
      statistic = Inf,
      p_value = 0,
      df = NA_real_,
      alpha = NA_real_,
      xi1 = NA_real_,
      xi2 = NA_real_)
    return(ret)
  } else if (sum(x > 0.5) == 1) {
    ## only one genotype, and data conform to it
    ret <- list(
      statistic = 0,
      p_value = 1,
      df = 0,
      alpha = NA_real_,
      xi1 = NA_real_,
      xi2 = NA_real_)
    return(ret)
  }

  if (pp && dr) {
    ret <- lrt_dr_pp_g4(x = x, g1 = g1, g2 = g2, drbound = drbound)
  } else if (pp && !dr) {
    ret <- lrt_ndr_pp_g4(x = x, g1 = g1, g2 = g2, alpha = alpha)
  } else if (!pp && dr) {
    ret <- lrt_dr_npp_g4(x = x, g1 = g1, g2 = g2, drbound = drbound, xi1 = xi1, xi2 = xi2)
  } else {
    ret <- lrt_ndr_npp_g4(x = x, g1 = g1, g2 = g2, alpha = alpha, xi1 = xi1, xi2 = xi2)
  }

  return(ret)
}

#' Checks to see if x is impossible given parental genotypes
#'
#' @param x The offspring genotype counts
#' @param g1 Parent 1's genotype
#' @param g2 Parent 2's genotype
#' @param df Is double reduction possible?
#'
#' @return TRUE if impossible, FALSE if possible
#'
#' @author David Gerard
#'
#' @examples
#' is_impossible(x = c(1, 2, 3, 0, 0), g1 = 0, g2 = 3, dr = TRUE)
#'
#'
#' @noRd
is_impossible <- function(x, g1, g2, dr = TRUE) {
  if (g1 == 0 && g2 == 0) {
    return(!all(x[2:5] == 0))
  } else if (g1 == 4 & g2 == 4) {
    return(!all(x[1:4] == 0))
  } else if (g1 == 0 && g2 == 4 || g1 == 4 && g2 == 0) {
    return(!all(x[c(1, 2, 4, 5)] == 0))
  } else if (dr) {
    if (g1 == 0 && g2 %in% c(1, 2, 3) || g1 %in% c(1, 2, 3) && g2 == 0) {
      return(!all(x[4:5] == 0))
    } else if (g1 == 4 && g2 %in% c(1, 2, 3) || g1 %in% c(1, 2, 3) && g2 == 4) {
      return(!all(x[1:2] == 0))
    }
  } else if (!dr) {
    if (g1 == 0 && g2 == 1 || g1 == 1 && g2 == 0) {
      return(!all(x[3:5] == 0))
    } else if (g1 == 0 && g2 == 2 || g1 == 2 && g2 == 0 || g1 == 1 && g2 == 1) {
      return(!all(x[4:5] == 0))
    } else if (g1 == 0 && g2 == 3 || g1 == 3 && g2 == 0) {
      return(!all(x[c(1, 4, 5)] == 0))
    } else if (g1 == 1 && g2 == 4 || g1 == 4 && g2 == 1) {
      return(!all(x[c(1, 2, 5)] == 0))
    } else if (g1 == 2 && g2 == 4 || g1 == 4 && g2 == 2 || g1 == 3 && g2 == 3) {
      return(!all(x[1:2] == 0))
    } else if (g1 == 3 && g2 == 4 || g1 == 4 && g2 == 3) {
      return(!all(x[1:3] == 0))
    } else if (g1 == 1 && g2 == 2 || g1 == 2 && g2 == 1) {
      return(!all(x[5] == 0))
    } else if (g1 == 1 && g2 == 3 || g1 == 3 && g2 == 1) {
      return(!all(x[c(1, 5)] == 0))
    } else if (g1 == 2 && g2 == 3 || g1 == 3 && g2 == 2) {
      return(!all(x[1] == 0))
    }
  }

  return(FALSE)
}

#' Returns number of theoretical zeros, to adjust df in lrt when genotypes are known
#'
#' @inheritParams is_impossible
#'
#' @noRd
nzeros_old <- function(g1, g2, dr = TRUE) {
  if (g1 == 0 && g2 == 0) {
    return(4)
  } else if (g1 == 4 & g2 == 4) {
    return(4)
  } else if (g1 == 0 && g2 == 4 || g1 == 4 && g2 == 0) {
    return(4)
  } else if (dr) {
    if (g1 == 0 && g2 %in% c(1, 2, 3) || g1 %in% c(1, 2, 3) && g2 == 0) {
      return(2)
    } else if (g1 == 4 && g2 %in% c(1, 2, 3) || g1 %in% c(1, 2, 3) && g2 == 4) {
      return(2)
    }
  } else if (!dr) {
    if (g1 == 0 && g2 == 1 || g1 == 1 && g2 == 0) {
      return(3)
    } else if (g1 == 0 && g2 == 2 || g1 == 2 && g2 == 0 || g1 == 1 && g2 == 1) {
      return(2)
    } else if (g1 == 0 && g2 == 3 || g1 == 3 && g2 == 0) {
      return(3)
    } else if (g1 == 1 && g2 == 4 || g1 == 4 && g2 == 1) {
      return(3)
    } else if (g1 == 2 && g2 == 4 || g1 == 4 && g2 == 2 || g1 == 3 && g2 == 3) {
      return(2)
    } else if (g1 == 3 && g2 == 4 || g1 == 4 && g2 == 3) {
      return(3)
    } else if (g1 == 1 && g2 == 2 || g1 == 2 && g2 == 1) {
      return(1)
    } else if (g1 == 1 && g2 == 3 || g1 == 3 && g2 == 1) {
      return(2)
    } else if (g1 == 2 && g2 == 3 || g1 == 3 && g2 == 2) {
      return(1)
    }
  }

  return(0)
}

#' Unlike nzeros_old, this one will check empirically for zeros
#'
#' @param g1 parent 1's genotype
#' @param g2 parent 2's genotype
#' @param alpha The double reduction rate (estimated or otherwise)
#' @param xi1 Parent 1's preferential pairing parameter (estimated or otherwise).
#' @param xi1 Parent 2's preferential pairing parameter (estimated or otherwise).
#'
#' @author David Gerard
#'
#' @noRd
nzeros <- function(g1, g2, alpha, xi1, xi2) {
  gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = g1, p2 = g2)
  TOL <- 1e-6
  return(sum(gf < TOL))
}


#' Tests if null parameters are estimated on boundary of parameterspace
#'
#' @inheritParams get_df
#'
#' @author David Gerard
#'
#' @return Either \code{TRUE} if the parameter is on the boundary, and
#' \code{FALSE} otherwise.
#'
#' @noRd
onbound <- function(g1, g2, alpha, xi1, xi2, dr, pp, drbound = 1/6, TOL = 1e-5) {
  ob <- TRUE
  if (dr && (g1 %in% c(1, 3)) && (g2 %in% c(1, 3))) {
    if (alpha > TOL && alpha < drbound - TOL) {
      ob <- FALSE
    } else {
      ## do nothing
    }
  } else if (dr && ((g1 %in% c(0, 4) && g2 %in% c(1, 3)) || (g1 %in% c(1, 3) && g2 %in% c(0, 4)))) {
    if (alpha > TOL && alpha < drbound - TOL) {
      ob <- FALSE
    } else {
      ## do nothing
    }
  } else if (dr && !pp && (g1 == 2 || g2 == 2)) {
    if (alpha > TOL && alpha < drbound - TOL) {
      ob <- FALSE
    } else {
      ## do nothing
    }
  } else if (pp && !dr && g1 == 2 && g2 != 2) {
    if (xi1 > TOL && xi1 < 1 - TOL) {
      ob <- FALSE
    } else {
      ## do nothing
    }
  } else if (pp && !dr && g1 != 2 && g2 == 2) {
    if (xi2 > TOL && xi2 < 1 - TOL) {
      ob <- FALSE
    } else {
      ## do nothing
    }
  } else if (pp && !dr && g1 == 2 && g2 == 2) {
    if ((xi1 < TOL || xi1 > 1 - TOL) && (xi2 < TOL || xi2 > 1 - TOL)) {
      ## do nothing
    } else {
      ob <- FALSE
    }
  } else if (pp && dr && g1 == 2 && g2 %in% c(0, 1, 3, 4)) {
    if (alpha < TOL && (xi1 < TOL || xi1 > 1 - TOL)) {
      ## do nothing
    } else if (alpha > drbound - TOL) {
      ## do nothing
    } else {
      ob <- FALSE
    }
  } else if (pp && dr && g1 %in% c(0, 1, 3, 4) && g2 == 2) {
    if (alpha < TOL && (xi2 < TOL || xi2 > 1 - TOL)) {
      ## do nothing
    } else if (alpha > drbound - TOL) {
      ## do nothing
    } else {
      ob <- FALSE
    }
  } else if (pp && dr && g1 == 2 && g2 == 2) {
    if (alpha < TOL && (xi1 < TOL || xi1 > 1 - TOL)) {
      ## do nothing
    } else if (alpha < TOL && (xi2 < TOL || xi2 > 1 - TOL)) {
      ## do nothing
    } else if (alpha > drbound - TOL) {
      ## do nothing
    } else {
      ob <- FALSE
    }
  } else {
    ## do nothing
    ## should only get here is !pp && !dr
  }

  return(ob)
}

#' Experimental calculation of degrees of freedom
#'
#' The complicated non-identifiability actually makes calculating the degrees
#' of freedom really hard. This function was created by experimentally running
#' through every possible scenario to see what is the empirical null
#' distribution. We are effectively using the strategy of Susko (2013), where
#' we estimate the degrees of freedom. The estimated number of degrees of freedom
#' under the alternative is the number of non-zero elements minus 1. The number
#' under the null is 1 if the parameter is not estimated on the boundary, and
#' zero otherwise.
#'
#' @param g1 Parent 1's genotype
#' @param g2 Parent 2's genotype
#' @param alpha The estimated double reduction rate
#' @param xi1 The estimated preferential pairing parameter for parent 1.
#' @param xi2 The estimated preferential pairing parmaeter for parent 2.
#' @param TOL tolerance to check boundary of parameter space.
#' @param dr Was double reduction being estimated?
#' @param pp Was preferential pairing being estimated?
#'
#' @author David Gerard
#'
#' @references Susko, E. (2013). Likelihood ratio tests with boundary constraints using data-dependent degrees of freedom. Biometrika, 100(4), 1019-1023.
#'
#' @noRd
get_df <- function(g1, g2, alpha, xi1, xi2, dr, pp, drbound = 1/6, TOL = 1e-5) {
  nz <- nzeros(g1 = g1, g2 = g2, alpha = alpha, xi1 = xi1, xi2 = xi2)

  ## Number of parameters under alternative
  df_alt <- 4 - nz ## Number of non-zeros (5 - nz) minuz 1, so 4 - nz

  df_null <- !onbound(
    g1 = g1,
    g2 = g2,
    alpha = alpha,
    xi1 = xi1,
    xi2 = xi2,
    dr = dr,
    pp = pp,
    drbound = drbound,
    TOL = TOL)

  df <- df_alt - df_null

  df <- max(df, 1)

  return(df)
}


#' Likelihood under three parameter model when genotypes are known
#'
#' This is under the three parameter model.
#'
#' @param x A vector of length 5. \code{x[i]} is the count of individuals with
#'     genotype \code{i-1}.
#' @param tau The probability of quadrivalent formation.
#' @param beta The probability of double reduction given quadrivalent formation.
#' @param gamma1 The probability of AA_aa pairing for parent 1.
#' @param gamma2 The probability of AA_aa pairing for parent 2.
#' @param g1 Parent 1's genotype.
#' @param g2 Parent 2's genotype.
#' @param log_p A logical. Should we return the log likelihood or not?
#' @param pen A tiny penalty to help with numerical stability
#'
#' @return The (log) likelihood.
#'
#' @examples
#' x <- c(1, 4, 5, 3, 1)
#' tau <- 0.5
#' beta <- 0.1
#' gamma1 <- 0.5
#' gamma2 <- 0.3
#' g1 <- 1
#' g2 <- 2
#' like_gknown_3(
#'   x = x,
#'   tau = tau,
#'   beta = beta,
#'   gamma1 = gamma1,
#'   gamma2 = gamma2,
#'   g1 = g1,
#'   g2 = g2)
#'
#' @author David Gerard
#'
#' @export
like_gknown_3 <- function(x, tau, beta, gamma1, gamma2, g1, g2, log_p = TRUE, pen = 0) {
  stopifnot(length(x) == 5)
  stopifnot(tau >= 0, tau <= 1,
            beta >= 0, beta <= 1,
            gamma1 >= 0, gamma1 <= 1,
            x >= 0,
            g1 >= 0, g1 <= 4,
            g2 >= 0, g2 <= 4)

  gf <- offspring_gf_3(
    tau = tau,
    beta = beta,
    gamma1 = gamma1,
    gamma2 = gamma2,
    p1 = g1,
    p2 = g2)

  gf <- gf + pen
  gf <- gf / sum(gf)

  return(stats::dmultinom(x = x, prob = gf, log = log_p))
}

#' Likelihood under three parameter model when genotypes are known
#'
#' This is under the two parameter model.
#'
#' @param x A vector of length 5. \code{x[i]} is the count of individuals with
#'     genotype \code{i-1}.
#' @param alpha The double reduction rate.
#' @param xi1 The preferential pairing parameter of parent 1.
#' @param xi2 The preferential pairing parameter of parent 2.
#' @param g1 Parent 1's genotype.
#' @param g2 Parent 2's genotype.
#' @param log_p A logical. Should we return the log likelihood or not?
#' @param pen A tiny penalty to help with numerical stability
#'
#' @return The (log) likelihood.
#'
#' @examples
#' x <- c(1, 4, 5, 3, 1)
#' alpha <- 0.01
#' xi1 <- 0.5
#' xi2 <- 0.3
#' g1 <- 1
#' g2 <- 2
#' like_gknown_2(
#'   x = x,
#'   alpha = alpha,
#'   xi1 = xi1,
#'   xi2 = xi2,
#'   g1 = g1,
#'   g2 = g2)
#'
#' @author David Gerard
#'
#' @export
like_gknown_2 <- function(x, alpha, xi1, xi2, g1, g2, log_p = TRUE, pen = 0) {
  stopifnot(length(x) == 5)
  stopifnot(alpha >= 0, alpha <= 1,
            xi1 >= 0, xi1 <= 1,
            xi2 >= 0, xi2 <= 1,
            x >= 0,
            g1 >= 0, g1 <= 4,
            g2 >= 0, g2 <= 4)
  gf <- offspring_gf_2(
    alpha = alpha,
    xi1 = xi1,
    xi2 = xi2,
    p1 = g1,
    p2 = g2)

  gf <- gf + pen
  gf <- gf / sum(gf)

  return(stats::dmultinom(x = x, prob = gf, log = log_p))
}

#' Likelihood ratio test assuming no double reduction and no preferential pairing
#'
#' This is when all genotypes are known.
#'
#' @inheritParams like_gknown_3
#' @param alpha The known rate of double reduction.
#' @param xi1 The known preferential pairing parameter of parent 1.
#' @param xi2 The known preferential pairing parameter of parent 2.
#'
#' @return A list of length three with the following elements
#' \describe{
#'   \item{\code{statistic}}{The likelihood ratio test statistic.}
#'   \item{\code{df}}{The degrees of freedom.}
#'   \item{\code{p_value}}{The p-value.}
#' }
#'
#' @author David Gerard
#'
#' @noRd
lrt_ndr_npp_g4 <- function(x, g1, g2, alpha = 0, xi1 = 1/3, xi2 = 1/3) {

  ## check impossible genotype ----
  TOL <- sqrt(.Machine$double.eps)
  if (is_impossible(x = x, g1 = g1, g2 = g2, dr = alpha > TOL)) {
    ret <- list(
      statistic = Inf,
      p_value = 0,
      df = NA_real_,
      alpha = NA_real_,
      xi1 = NA_real_,
      xi2 = NA_real_)
    return(ret)
  }


  ## Run test
  l1 <- stats::dmultinom(x = x, prob = x / sum(x), log = TRUE)
  l0 <- like_gknown_2(
    x = x,
    alpha = alpha,
    xi1 = xi1,
    xi2 = xi2,
    g1 = g1,
    g2 = g2,
    log_p = TRUE)

  llr <- -2 * (l0 - l1)

  ## Get degrees of freedom ----
  df <- get_df(g1 = g1, g2 = g2, alpha = alpha, xi1 = xi1, xi2 = xi2, dr = FALSE, pp = FALSE)

  p_value <- stats::pchisq(q = llr, df = df, lower.tail = FALSE, log.p = FALSE)

  ret <- list(
    statistic = llr,
    p_value = p_value,
    df = df,
    alpha = alpha,
    xi1 = xi1,
    xi2 = xi2)

  return(ret)
}

## Methods when dr and pp are not known ---------------------------------------

#' Objective when estimating both double reduction and preferential pairing
#'
#' @param par either just alpha (when no parent genotype is 2)
#'   or (tau, beta, gamma1) or (tau, beta, gamma2) or
#'   (tau, beta, gamma1, gamma2).
#' @param x offspring genotype counts
#' @param g1 first parent genotype
#' @param g2 second parent genotype
#' @param pen A tiny penalty to help with numerical stability
#'
#' @author David Gerard
#'
#' @noRd
obj_dr_pp <- function(par, x, g1, g2, pen = 1e-6) {
  par[par < 0] <- 0 ## to deal with -1e-18 etc
  par[par > 1] <- 1
  if (g1 != 2 && g2 != 2) {
    stopifnot(length(par) == 1)
    obj <- like_gknown_3(
      x = x,
      tau = 1,
      beta = par[[1]],
      gamma1 = 1/3,
      gamma2 = 1/3,
      g1 = g1,
      g2 = g2,
      log_p = TRUE,
      pen = pen)
  } else if (g1 == 2 && g2 != 2){
    stopifnot(length(par) == 3)
    obj <- like_gknown_3(
      x = x,
      tau = par[[1]],
      beta = par[[2]],
      gamma1 = par[[3]],
      gamma2 = 1/3,
      g1 = g1,
      g2 = g2,
      log_p = TRUE,
      pen = pen)
  } else if (g1 != 2 && g2 == 2){
    stopifnot(length(par) == 3)
    obj <- like_gknown_3(
      x = x,
      tau = par[[1]],
      beta = par[[2]],
      gamma1 = 1/3,
      gamma2 = par[[3]],
      g1 = g1,
      g2 = g2,
      log_p = TRUE,
      pen = pen)
  } else {
    stopifnot(length(par) == 4)
    obj <- like_gknown_3(
      x = x,
      tau = par[[1]],
      beta = par[[2]],
      gamma1 = par[[3]],
      gamma2 = par[[4]],
      g1 = g1,
      g2 = g2,
      log_p = TRUE,
      pen = pen)
  }
  return(obj)
}

#' Initial parameter values
#'
#' @param g1 parent 1 genotype
#' @param g2 parent 2 genotype.
#' @param drbound The maximum double reduction rate
#' @param type How should we initialize? "random" or "half"?
#' @param fudge How much to add to the lower bound and subtract from
#'     the upper bound.
#'
#' @return A list of length 3. The first element is the initial values,
#'     the second is the lower bound of parameters, the third element
#'     is the upper bound of parameters.
#'
#' @author David Gerard
#'
#' @noRd
lrt_init <- function(g1, g2, drbound = 1/6, type = c("random", "half"), fudge = 0) {
  type <- match.arg(type)
  if (type == "random") {
    mult <- stats::runif(4)
  } else if (type == "half") {
    mult <- c(0.5, 0.5, 1/3, 1/3)
  }

  if (g1 != 2 && g2 != 2) {
    out <- list(
      par = c(alpha = drbound * mult[[1]]),
      lower = fudge,
      upper = drbound)
  } else if (g1 == 2 & g2 != 2) {
    out <- list(
      par = c(tau = mult[[1]], beta = drbound * mult[[2]], gamma1 = mult[[3]]),
      lower = rep(fudge, length.out = 3),
      upper = c(1 - fudge, drbound, 1 - fudge))
  } else if (g1 != 2 & g2 == 2) {
    out <- list(
      par = c(tau = mult[[1]], beta = drbound * mult[[2]], gamma2 = mult[[3]]),
      lower = rep(fudge, length.out = 3),
      upper = c(1 - fudge, drbound, 1 - fudge))
  } else {
    out <- list(
      par = c(tau = mult[[1]], beta = drbound * mult[[2]], gamma1 = mult[[3]], gamma2 = mult[[4]]),
      lower = rep(fudge, length.out = 4),
      upper = c(1 - fudge, drbound, 1 - fudge, 1 - fudge))
  }
  return(out)
}

#' LRT when both double reduction and preferential pairing are not known.
#'
#' All genotypes are known.
#'
#' @inheritParams like_gknown_3
#' @param drbound the upper bound on the double reduction rate.
#' @param ntry The number of times to run the gradient ascent.
#'
#' @author David Gerard
#'
#' @examples
#' set.seed(1000)
#' alpha <- 1/12
#' xi1 <- 1/3
#' xi2 <- 1/3
#' n <- 1000
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 0, p2 = 0)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 0, g2 = 0)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 1, p2 = 0)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 1, g2 = 0)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 1, p2 = 1)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 1, g2 = 1)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 1, p2 = 2)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 1, g2 = 2)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 1, p2 = 3)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 1, g2 = 3)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 1, p2 = 4)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 1, g2 = 4)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 2, p2 = 0)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 2, g2 = 0)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 2, p2 = 1)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 2, g2 = 1)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 2, p2 = 2)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 2, g2 = 2)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 2, p2 = 3)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 2, g2 = 3)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 2, p2 = 4)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 2, g2 = 4)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 3, p2 = 0)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 3, g2 = 0)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 3, p2 = 1)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 3, g2 = 1)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 3, p2 = 2)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 3, g2 = 2)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 3, p2 = 3)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 3, g2 = 3)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 3, p2 = 4)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 3, g2 = 4)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 4, p2 = 4)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 4, g2 = 4)
#'
#' @noRd
lrt_dr_pp_g4 <- function(x, g1, g2, drbound = 1/6, ntry = 5) {

  if (g1 %in% c(0, 4) && g2 %in% c(0, 4)) {
    ret <- lrt_ndr_npp_g4(x = x, g1 = g1, g2 = g2)
    return(ret)
  }

  ## Deal with impossible values under null ----
  if (is_impossible(x = x, g1 = g1, g2 = g2, dr = TRUE)) {
    ret <- list(
      statistic = Inf,
      p_value = 0,
      df = NA_real_,
      alpha = NA_real_,
      xi1 = NA_real_,
      xi2 = NA_real_)
    return(ret)
  }

  ## max likelihood under alt
  l1 <- stats::dmultinom(x = x, prob = x / sum(x), log = TRUE)

  ## Find MLE under null
  l0 <- -Inf
  bout <- NULL
  for (i in seq_len(ntry)) {
    params <- lrt_init(g1 = g1, g2 = g2, drbound = drbound, type = "random")
    method <- ifelse(length(params$par) == 1, "Brent", "L-BFGS-B")
    oout <- stats::optim(
      par = params$par,
      fn = obj_dr_pp,
      method = method,
      lower = params$lower,
      upper = params$upper,
      control = list(fnscale = -1),
      x = x,
      g1 = g1,
      g2 = g2)

    if (oout$value > l0) {
      l0 <- oout$value
      bout <- oout
    }
  }

  ## Get df and calculate p-value
  if (g1 != 2 & g2 != 2) {
    alpha <- bout$par[[1]]
    xi1 <- NA_real_
    xi2 <- NA_real_
  } else if (g1 == 2 & g2 != 2) {
    tau <- bout$par[[1]]
    beta <- bout$par[[2]]
    gamma1 <- bout$par[[3]]
    two <- three_to_two(tau = tau, beta = beta, gamma = gamma1)
    alpha <- two[[1]]
    xi1 <- two[[2]]
    xi2 <- NA_real_
  } else if (g1 != 2 & g2 == 2) {
    tau <- bout$par[[1]]
    beta <- bout$par[[2]]
    gamma2 <- bout$par[[3]]
    two <- three_to_two(tau = tau, beta = beta, gamma = gamma2)
    alpha <- two[[1]]
    xi1 <- NA_real_
    xi2 <- two[[2]]
  } else if (g1 == 2 & g2 == 2) {
    tau <- bout$par[[1]]
    beta <- bout$par[[2]]
    gamma1 <- bout$par[[3]]
    gamma2 <- bout$par[[4]]
    two <- three_to_two(tau = tau, beta = beta, gamma = gamma1)
    alpha <- two[[1]]
    xi1 <- two[[2]]
    two <- three_to_two(tau = tau, beta = beta, gamma = gamma2)
    stopifnot(alpha == two[[1]])
    xi2 <- two[[2]]
  }

  ## Get degrees of freedom ----
  df <- get_df(g1 = g1, g2 = g2, alpha = alpha, xi1 = xi1, xi2 = xi2, dr = TRUE, pp = TRUE, drbound = drbound)

  llr <- -2 * (l0 - l1)

  p_value <- stats::pchisq(q = llr, df = df, lower.tail = FALSE)

  ret <-  list(
    statistic = llr,
    p_value = p_value,
    df = df,
    alpha = alpha,
    xi1 = xi1,
    xi2 = xi2)

  return(ret)
}

## No dr and estimate pp ------------------------------------------------------

#' Objective for lrt_ndr_pp_g4
#'
#' @inheritParams obj_dr_pp
#' @param par if g1 != 2 or g2 != 2, then should be xi1 (or xi2). If
#'     both g1 == 2 and g2 == 2, then first element is xi1 and second
#'     element is xi2.
#' @param alpha The known rate of double reduction.
#' @param pen A tiny penalty to help with numerical stability.
#'
#' @author David Gerard
#'
#' @noRd
obj_ndr_pp <- function(par, x, g1, g2, alpha = 0, pen = 1e-6) {
  par[par < 0] <- 0 ## to deal with -1e-18 etc
  par[par > 1] <- 1
  if (g1 != 2 && g2 != 2) {
    stop("obj_ndr_pp: have to have g1 == 2 or g2 == 2")
  } else if (g1 == 2 && g2 != 2) {
    stopifnot(length(par) == 1)
    obj <- like_gknown_2(
      x = x,
      alpha = alpha,
      xi1 = par[[1]],
      xi2 = 1/3,
      g1 = g1,
      g2 = g2,
      log_p = TRUE,
      pen = pen)
  } else if (g1 != 2 && g2 == 2) {
    stopifnot(length(par) == 1)
    obj <- like_gknown_2(
      x = x,
      alpha = alpha,
      xi1 = 1/3,
      xi2 = par[[1]],
      g1 = g1,
      g2 = g2,
      log_p = TRUE,
      pen = pen)
  } else if (g1 == 2 && g2 == 2) {
    stopifnot(length(par) == 2)
    obj <- like_gknown_2(
      x = x,
      alpha = alpha,
      xi1 = par[[1]],
      xi2 = par[[2]],
      g1 = g1,
      g2 = g2,
      log_p = TRUE,
      pen = pen)
  }
  return(obj)
}

#' LRT when double reduction is known, but not preferential pairing.
#'
#' All genotypes are known.
#'
#' @inherit lrt_dr_pp_g4
#' @param alpha The known rate of double reduction.
#' @param fudge How much to add to lower bound or subtract from upper bound.
#'
#' @author David Gerard
#'
#' @examples
#' gf <- offspring_gf_2(alpha = 0, xi1 = 1/3, xi2 = 1/3, p1 = 0, p2 = 2)
#' x <- offspring_geno(gf = gf, n = 100)
#' lrt_ndr_pp_g4(x = x, g1 = 0, g2 = 2)
#'
#' gf <- offspring_gf_2(alpha = 0, xi1 = 1/3, xi2 = 1/3, p1 = 1, p2 = 2)
#' x <- offspring_geno(gf = gf, n = 100)
#' lrt_ndr_pp_g4(x = x, g1 = 1, g2 = 2)
#'
#' gf <- offspring_gf_2(alpha = 0, xi1 = 1/3, xi2 = 1/3, p1 = 2, p2 = 2)
#' x <- offspring_geno(gf = gf, n = 100)
#' lrt_ndr_pp_g4(x = x, g1 = 2, g2 = 2)
#'
#' gf <- offspring_gf_2(alpha = 0, xi1 = 1/3, xi2 = 1/3, p1 = 3, p2 = 2)
#' x <- offspring_geno(gf = gf, n = 100)
#' lrt_ndr_pp_g4(x = x, g1 = 3, g2 = 2)
#'
#' gf <- offspring_gf_2(alpha = 0, xi1 = 1/3, xi2 = 1/3, p1 = 4, p2 = 2)
#' x <- offspring_geno(gf = gf, n = 100)
#' lrt_ndr_pp_g4(x = x, g1 = 4, g2 = 2)
#'
#' gf <- offspring_gf_2(alpha = 0, xi1 = 1/3, xi2 = 1/3, p1 = 2, p2 = 0)
#' x <- offspring_geno(gf = gf, n = 100)
#' lrt_ndr_pp_g4(x = x, g1 = 2, g2 = 0)
#'
#' gf <- offspring_gf_2(alpha = 0, xi1 = 1/3, xi2 = 1/3, p1 = 2, p2 = 1)
#' x <- offspring_geno(gf = gf, n = 100)
#' lrt_ndr_pp_g4(x = x, g1 = 2, g2 = 1)
#'
#' gf <- offspring_gf_2(alpha = 0, xi1 = 1/3, xi2 = 1/3, p1 = 2, p2 = 3)
#' x <- offspring_geno(gf = gf, n = 100)
#' lrt_ndr_pp_g4(x = x, g1 = 2, g2 = 3)
#'
#' gf <- offspring_gf_2(alpha = 0, xi1 = 1/3, xi2 = 1/3, p1 = 2, p2 = 4)
#' x <- offspring_geno(gf = gf, n = 100)
#' lrt_ndr_pp_g4(x = x, g1 = 2, g2 = 4)
#'
#' @noRd
lrt_ndr_pp_g4 <- function(x, g1, g2, alpha = 0, ntry = 5, fudge = 0) {
  if (g1 != 2 && g2 != 2) {
    return(lrt_ndr_npp_g4(x = x, g1 = g1, g2 = g2, alpha = alpha))
  }

  ## look at impossible genotypes under null ----
  TOL <- sqrt(.Machine$double.eps)
  if (is_impossible(x = x, g1 = g1, g2 = g2, dr = alpha > TOL)) {
    ret <- list(
      statistic = Inf,
      p_value = 0,
      df = NA_real_,
      alpha = NA_real_,
      xi1 = NA_real_,
      xi2 = NA_real_)
    return(ret)
  }

  ## run test ----

  ## max likelihood under alt
  l1 <- stats::dmultinom(x = x, prob = x / sum(x), log = TRUE)

  ## max likelihood under null
  if (g1 == 2 && g2 != 2) {
    oout <- stats::optim(
      par = 1/3,
      fn = obj_ndr_pp,
      lower = fudge,
      upper = 1 - fudge,
      method = "Brent",
      control = list(fnscale = -1),
      x = x,
      alpha = alpha,
      g1 = g1,
      g2 = g2)
    xi1 <- oout$par[[1]]
    xi2 <- NA_real_
  } else if (g1 != 2 && g2 == 2) {
    oout <- stats::optim(
      par = 1/3,
      fn = obj_ndr_pp,
      lower = fudge,
      upper = 1 - fudge,
      method = "Brent",
      control = list(fnscale = -1),
      x = x,
      alpha = alpha,
      g1 = g1,
      g2 = g2)
    xi1 <- NA_real_
    xi2 <- oout$par[[1]]
  } else if (g1 == 2 && g2 == 2) {
    oout <- list()
    oout$value <- -Inf
    for (i in seq_len(ntry)) {
      par_init <- stats::runif(n = 2, min = 0, max = 1)
      bout <- stats::optim(
        par = par_init,
        fn = obj_ndr_pp,
        lower = c(fudge, fudge),
        upper = c(1 - fudge, 1 - fudge),
        method = "L-BFGS-B",
        control = list(fnscale = -1),
        x = x,
        alpha = alpha,
        g1 = g1,
        g2 = g2)
      if (bout$value > oout$value) {
        oout <- bout
      }
    }
    xi1 <- oout$par[[1]]
    xi2 <- oout$par[[2]]
  }
  l0 <- oout$value

  ## Get degrees of freedom ----
  df <- get_df(g1 = g1, g2 = g2, alpha = alpha, xi1 = xi1, xi2 = xi2, dr = FALSE, pp = TRUE)

  llr <- -2 * (l0 - l1)

  p_value <- stats::pchisq(q = llr, df = df, lower.tail = FALSE)

  ret <-  list(
    statistic = llr,
    p_value = p_value,
    df = df,
    alpha = alpha,
    xi1 = xi1,
    xi2 = xi2)

  return(ret)
}

## Estimate DR, fix pp -------------------------------------------------------


#' Wrapper for like_gknnown_2 with rounding at 0 and 1.
#'
#' @inheritParams like_gknown_2
#' @param par Alpha, the double reduction rate
#'
#' @author David Gerard
#'
#' @noRd
obj_dr_npp <- function(par, x, xi1, xi2, g1, g2, pen) {
  par[par < 0] <- 0 ## deal with -1e-18 etc
  par[par > 1] <- 1
  ll <- like_gknown_2(x = x, alpha = par, xi1 = xi1, xi2 = xi2, g1 = g1, g2 = g2, log_p = TRUE, pen = pen)
  return(ll)
}

#' LRT when double reduction is not known, but no preferential pairing.
#'
#' All genotypes are known.
#'
#' @inherit lrt_dr_pp_g4
#' @param xi1 The known preferential pairingn parameter of parent 1.
#' @param xi2 The known preferential pairingn parameter of parent 2.
#' @param fudge How much to add to lower bound or subtract from upper bound.
#'
#' @author David Gerard
#'
#' @examples
#'
#' gf <- offspring_gf_2(alpha = 1/12, xi1 = 1/3, xi2 = 1/3, p1 = 0, p2 = 2)
#' x <- offspring_geno(gf = gf, n = 100)
#' lrt_dr_npp_g4(x = x, g1 = 0, g2 = 2)
#'
#' gf <- offspring_gf_2(alpha = 1/12, xi1 = 1/3, xi2 = 1/3, p1 = 1, p2 = 2)
#' x <- offspring_geno(gf = gf, n = 100)
#' lrt_dr_npp_g4(x = x, g1 = 1, g2 = 2)
#'
#' gf <- offspring_gf_2(alpha = 1/12, xi1 = 1/3, xi2 = 1/3, p1 = 2, p2 = 2)
#' x <- offspring_geno(gf = gf, n = 100)
#' lrt_dr_npp_g4(x = x, g1 = 2, g2 = 2)
#'
#' gf <- offspring_gf_2(alpha = 1/12, xi1 = 1/3, xi2 = 1/3, p1 = 3, p2 = 2)
#' x <- offspring_geno(gf = gf, n = 100)
#' lrt_dr_npp_g4(x = x, g1 = 3, g2 = 2)
#'
#' gf <- offspring_gf_2(alpha = 1/12, xi1 = 1/3, xi2 = 1/3, p1 = 4, p2 = 2)
#' x <- offspring_geno(gf = gf, n = 100)
#' lrt_dr_npp_g4(x = x, g1 = 4, g2 = 2)
#'
#' @noRd
lrt_dr_npp_g4 <- function(x, g1, g2, drbound = 1/6, xi1 = 1/3, xi2 = 1/3, fudge = 0) {
  ## Same as in lrt_dr_pp_g4 when no 2's
  if (g1 != 2 & g2 != 2) {
    return(lrt_dr_pp_g4(x = x, g1 = g1, g2 = g2, drbound = drbound))
  }

  ## Deal with impossible genotypes under null
  if (is_impossible(x = x, g1 = g1, g2 = g2, dr = TRUE)) {
    ret <- list(
      statistic = Inf,
      p_value = 0,
      df = NA_real_,
      alpha = NA_real_,
      xi1 = NA_real_,
      xi2 = NA_real_)
    return(ret)
  }

  ## max likelihood under alt
  l1 <- stats::dmultinom(x = x, prob = x / sum(x), log = TRUE)

  ## max likelihood under null
  oout <- stats::optim(par = drbound / 2,
                       fn = obj_dr_npp,
                       method = "Brent",
                       lower = fudge,
                       upper = drbound,
                       control = list(fnscale = -1),
                       x = x,
                       xi1 = xi1,
                       xi2 = xi2,
                       g1 = g1,
                       g2 = g2,
                       pen = 1e-6)
  l0 <- oout$value
  alpha <- oout$par[[1]]

  ## Get degrees of freedom ----
  df <- get_df(g1 = g1, g2 = g2, alpha = alpha, xi1 = xi1, xi2 = xi2, dr = TRUE, pp = FALSE, drbound = drbound)

  llr <- -2 * (l0 - l1)
  p_value <- stats::pchisq(q = llr, df = df, lower.tail = FALSE)

  ret <- list(
    statistic = llr,
    p_value = p_value,
    df = df,
    alpha = alpha,
    xi1 = xi1,
    xi2 = xi2)

  return(ret)
}
