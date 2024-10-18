## Offspring genotype likelihoods, parent genotypes

#' Likelihood ratio test using genotype likelihoods.
#'
#' This will run a likelihood ratio test using the genotypes of an F1 population
#' of tetraploids for the null of Mendelian segregation (accounting for double
#' reduction and preferential pairing) against the alternative of
#' segregation distortion. This is when genotype uncertainty is accounted
#' for through genotype likelihoods.
#'
#' @inheritSection lrt_men_g4 Unidentified parameters
#'
#' @inherit lrt_men_g4 return
#'
#' @param gl The genotype log-likelihoods. The rows index the individuals
#'    and the columns index the genotypes.
#' @param g1 Either parent 1's genotype, or parent 1's genotype log-likelihoods.
#' @param g2 Either parent 2's genotype, or parent 2's genotype log-likelihoods.
#' @param drbound The upper bound on the double reduction rate.
#' @param pp Is (partial) preferential pairing possible (\code{TRUE}) or
#'    not (\code{FALSE})?
#' @param dr Is double reduction possible (\code{TRUE}) or
#'    not (\code{FALSE})?
#' @param alpha If \code{dr = FALSE}, this is the known rate of double
#'     reduction.
#' @param xi1 If \code{pp = FALSE}, this is the known preferential pairing
#'     parameter of parent 1.
#' @param xi2 If \code{pp = FALSE}, this is the known preferential pairing
#'     parameter of parent 2.
#'
#' @author David Gerard
#'
#' @examples
#' ## null simulation
#' set.seed(1)
#' g1 <- 2
#' g2 <- 2
#' gl <- simf1gl(n = 25, g1 = g1, g2 = g2, alpha = 1/12, xi2 = 1/2)
#'
#' ## LRT when parent genotypes are known.
#' lrt_men_gl4(gl = gl, g1 = g1, g2 = g2)
#'
#' ## LRT when parent genotypes are not known
#' lrt_men_gl4(gl = gl)
#'
#' ## Alternative simulation
#' gl <- simgl(nvec = rep(5, 5))
#' lrt_men_gl4(gl = gl, g1 = g1, g2 = g2)
#'
#' @export
lrt_men_gl4 <- function(
    gl,
    g1 = NULL,
    g2 = NULL,
    drbound = 1/6,
    pp = TRUE,
    dr = TRUE,
    alpha = 0,
    xi1 = 1/3,
    xi2 = 1/3) {
  stopifnot(ncol(gl) == 5,
            length(drbound) == 1,
            drbound <= 1,
            drbound > 1e-6,
            alpha >= 0,
            alpha <= 1,
            xi1 >= 0,
            xi1 <= 1,
            xi2 >= 0,
            xi2 <= 1,
            length(pp) == 1,
            is.logical(pp),
            length(dr) == 1,
            is.logical(dr),
            length(alpha) == 1,
            length(xi1) == 1,
            length(xi2) == 1)

  ## if parent geno is NULL, then use unknown genotypes method
  if (is.null(g1) && is.null(g2)) {
    g1 <- rep(-log(5), 5)
    g2 <- rep(-log(5), 5)
  } else if (is.null(g1) && !is.null(g2)) {
    g1 <- rep(-log(5), 5)
    if (length(g2) == 1) {
      g2temp <- rep(-Inf, 5)
      g2temp[[g2 + 1]] <- 0
      g2 <- g2temp
    }
  } else if (!is.null(g1) && is.null(g2)) {
    g2 <- rep(-log(5), 5)
    if (length(g1) == 1) {
      g1temp <- rep(-Inf, 5)
      g1temp[[g1 + 1]] <- 0
      g1 <- g1temp
    }
  }

  ## Check parent genotypes
  if (length(g1) == 1 && length(g2) == 1) {
    pknown <- TRUE
    g1 <- round(g1)
    g2 <- round(g2)
    stopifnot(g1 >= 0, g1 <= 4, g2 >= 0, g2 <= 4)
  } else if (length(g1) == 5 && length(g2) == 5) {
    pknown <- FALSE
    ## normalize to sum to 1
    g1 <- g1 - log_sum_exp(g1)
    g2 <- g2 - log_sum_exp(g2)
  } else {
    stop("g1 and g2 should either both be length 1 (known genotypes)\nor both be length 5 (genotype log-likelihoods).")
  }

  ## run test
  if (pp && dr && pknown) {
    ret <- lrt_dr_pp_glpknown4(gl = gl, g1 = g1, g2 = g2, drbound = drbound)
  } else if (pp && !dr && pknown) {
    ret <- lrt_ndr_pp_glpknown4(gl = gl, g1 = g1, g2 = g2, alpha = alpha)
  } else if (!pp && dr && pknown) {
    ret <- lrt_dr_npp_glpknown4(gl = gl, g1 = g1, g2 = g2, drbound = drbound, xi1 = xi1, xi2 = xi2)
  } else if (!pp && !dr && pknown) {
    ret <- lrt_ndr_npp_glpknown4(gl = gl, g1 = g1, g2 = g2, alpha = alpha, xi1 = xi1, xi2 = xi2)
  } else {
    ret <- lrt_gl4(
      gl = gl,
      p1_gl = g1,
      p2_gl = g2,
      drbound = drbound,
      dr = dr,
      pp = pp,
      alpha = alpha,
      xi1 = xi1,
      xi2 = xi2)
  }

  return(ret)
}

#' Likelihood under three parameter model when using offspring genotypes
#' likelihoods but parent genotypes are known.
#'
#' This is under the three parameter model.
#'
#' @inheritParams like_gknown_3
#' @param gl The matrix of genotype likelihoods of the offspring. Rows index
#'     The individuals, columns index the genotypes.
#'
#' @author David Gerard
#'
#' @return The (log) likelihood of the three parameter model when using
#'     genotype likelihoods.
#'
#' @examples
#' g1 <- 1
#' g2 <- 0
#' gl <- simf1gl(
#'   n = 25,
#'   g1 = g1,
#'   g2 = g2,
#'   rd = 10,
#'   alpha = 0,
#'   xi1 = 1/3,
#'   xi2 = 1/3)
#' like_glpknown_3(
#'   gl = gl,
#'   tau = 1/2,
#'   beta = 1/12,
#'   gamma1 = 1/3,
#'   gamma2 = 1/3,
#'   g1 = g1,
#'   g2 = g2,
#'   log_p = TRUE)
#'
#' @export
like_glpknown_3 <- function(gl, tau, beta, gamma1, gamma2, g1, g2, log_p = TRUE) {
  stopifnot(ncol(gl) == 5)
  stopifnot(tau >= 0, tau <= 1,
            beta >= 0, beta <= 1,
            gamma1 >= 0, gamma1 <= 1,
            g1 >= 0, g1 <= 4,
            g2 >= 0, g2 <= 4)

  gf <- offspring_gf_3(
    tau = tau,
    beta = beta,
    gamma1 = gamma1,
    gamma2 = gamma2,
    p1 = g1,
    p2 = g2)

  return(llike_li(B = gl, lpivec = log(gf)))
}

#' Likelihood under three parameter model when using offspring genotypes
#' likelihoods but parent genotypes are known.
#'
#' This is under the two parameter model.
#'
#' @inheritParams like_gknown_2
#' @param gl The matrix of genotype likelihoods of the offspring. Rows index
#'     The individuals, columns index the genotypes.
#'
#' @author David Gerard
#'
#' @return The (log) likelihood of the two parameter model when using
#'     genotype likelihoods.
#'
#' @examples
#' g1 <- 1
#' g2 <- 0
#' gl <- simf1gl(
#'   n = 25,
#'   g1 = g1,
#'   g2 = g2,
#'   rd = 10,
#'   alpha = 0,
#'   xi1 = 1/3,
#'   xi2 = 1/3)
#' like_glpknown_2(
#'   gl = gl,
#'   alpha = 0.01,
#'   xi1 = 0.5,
#'   xi2 = 0.3,
#'   g1 = g1,
#'   g2 = g2,
#'   log_p = TRUE)
#'
#' @export
like_glpknown_2 <- function(gl, alpha, xi1, xi2, g1, g2, log_p = TRUE) {
  stopifnot(ncol(gl) == 5)
  stopifnot(alpha >= 0, alpha <= 1,
            xi1 >= 0, xi1 <= 1,
            xi2 >= 0, xi2 <= 1,
            g1 >= 0, g1 <= 4,
            g2 >= 0, g2 <= 4)

  gf <- offspring_gf_2(
    alpha = alpha,
    xi1 = xi1,
    xi2 = xi2,
    p1 = g1,
    p2 = g2)

  return(llike_li(B = gl, lpivec = log(gf)))
}


#' Likelihood ratio test assuming no double reduction and no preferential pairing
#'
#' This is when offspring genotypes are not known
#'
#' @inheritParams like_glpknown_3
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
lrt_ndr_npp_glpknown4 <- function(gl, g1, g2, alpha = 0, xi1 = 1/3, xi2 = 1/3) {
  ## MLE under alternative
  log_qhat1 <- c(em_li(B = gl))
  l1 <- llike_li(B = gl, lpivec = log_qhat1)

  ## null likelihood
  qhat2 <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = g1, p2 = g2)
  l0 <- llike_li(B = gl, lpivec = log(qhat2))

  ## calculate test statistic
  llr <- -2 * (l0 - l1)

  ## calculate degrees of freedom
  TOL <- 1e-3
  df_alt <- sum(offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = g1, p2 = g2) > TOL | log_qhat1 > log(TOL)) - 1
  df_null <- 0 ## fixed
  df <- df_alt - df_null
  df <- max(df, 1)

  ## calculate p-value
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

## LRT when dr and pp are estimated --------------------------------------

#' Objective for lrt_dr_pp_glpknown4
#'
#' @param par either just alpha (when no parent genotype is 2)
#'   or (tau, beta, gamma1) or (tau, beta, gamma2) or
#'   (tau, beta, gamma1, gamma2).
#' @param x offspring genotype counts
#' @param g1 first parent genotype
#' @param g2 second parent genotype
#'
#' @return log likelihood
#'
#' @author David Gerard
#' @noRd
obj_dr_pp_gl <- function(par, gl, g1, g2) {
  par[par < 0] <- 0 ## deal with -1e-18 etc
  par[par > 1] <- 1
  if (g1 != 2 && g2 != 2) {
    stopifnot(length(par) == 1)
    obj <- like_glpknown_3(
      gl = gl,
      tau = 1,
      beta = par[[1]],
      gamma1 = 1/3,
      gamma2 = 1/3,
      g1 = g1,
      g2 = g2)
  } else if (g1 == 2 && g2 != 2){
    stopifnot(length(par) == 3)
    obj <- like_glpknown_3(
      gl = gl,
      tau = par[[1]],
      beta = par[[2]],
      gamma1 = par[[3]],
      gamma2 = 1/3,
      g1 = g1,
      g2 = g2)
  } else if (g1 != 2 && g2 == 2){
    stopifnot(length(par) == 3)
    obj <- like_glpknown_3(
      gl = gl,
      tau = par[[1]],
      beta = par[[2]],
      gamma1 = 1/3,
      gamma2 = par[[3]],
      g1 = g1,
      g2 = g2)
  } else {
    stopifnot(length(par) == 4)
    obj <- like_glpknown_3(
      gl = gl,
      tau = par[[1]],
      beta = par[[2]],
      gamma1 = par[[3]],
      gamma2 = par[[4]],
      g1 = g1,
      g2 = g2)
  }
  return(obj)
}

#' LRT when both double reduction and preferential pairing are not known.
#'
#' Offspring use genotype likelihoods, parent genotypes are known.
#'
#' @inheritParams like_glpknown_3
#' @param drbound Maximum value of double reduction.
#' @param ntry The number of times to run the gradient ascent.
#'
#' @author David Gerard
#'
#' @examples
#' ## null sim
#' set.seed(10)
#' g1 <- 2
#' g2 <- 2
#' gl <- simf1gl(n = 25, g1 = g1, g2 = g2, alpha = 1/12, xi2 = 1/4)
#' lrt_dr_pp_glpknown4(gl = gl, g1 = g1, g2 = g2)
#'
#' ## Alt sim
#' gl <- simgl(nvec = rep(5, 5))
#' lrt_dr_pp_glpknown4(gl = gl, g1 = g1, g2 = g2)
#'
#' @noRd
lrt_dr_pp_glpknown4 <- function(gl, g1, g2, drbound = 1/6, ntry = 5) {
  ## MLE under alternative
  log_qhat1 <- c(em_li(B = gl))
  l1 <- llike_li(B = gl, lpivec = log_qhat1)

  ## MLE under Null
  l0 <- -Inf
  bout <- NULL
  for (i in seq_len(ntry)) {
    params <- lrt_init(g1 = g1, g2 = g2, drbound = drbound, type = "random")
    method <- ifelse(length(params$par) == 1, "Brent", "L-BFGS-B")
    oout <- stats::optim(
      par = params$par,
      fn = obj_dr_pp_gl,
      method = method,
      lower = params$lower,
      upper = params$upper,
      control = list(fnscale = -1),
      gl = gl,
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

  ## calculate degrees of freedom
  TOL <- 1e-3
  df_alt <- sum(offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = g1, p2 = g2) > TOL | log_qhat1 > log(TOL)) - 1
  df_null <- !onbound(
    g1 = g1,
    g2 = g2,
    alpha = alpha,
    xi1 = xi1,
    xi2 = xi2,
    dr = TRUE,
    pp = TRUE,
    drbound = drbound,
    TOL = TOL)
  df <- df_alt - df_null
  df <- max(df, 1)

  ## calculate LRT stat
  llr <- -2 * (l0 - l1)

  ## calculate p-value
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


## No pp ----------------------------------------------------------------------

#' Wrapper for like_glpknown_2 with rounding at 0 and -1
#'
#' @inheritParams like_glpknown_2
#' @param par The double reduction rate alpha
#'
#' @author David Gerard
#'
#' @noRd
obj_dr_npp_gl <- function(par, gl, xi1, xi2, g1, g2) {
  par[par < 0] <- 0
  par[par > 1] <- 1
  ll <- like_glpknown_2(gl = gl, alpha = par, xi1 = xi1, xi2 = xi2, g1 = g1, g2 = g2, log_p = TRUE)
  return(ll)
}

#' Likelihood ratio test, gl for offspring, known for parents, dr, no pp.
#'
#' @inheritParams lrt_dr_pp_glpknown4
#' @param xi1 The known preferential pairing parameter of parent 1.
#' @param xi2 The known preferential pairing parameter of parent 2.
#' @param fudge How much to add to the lower bound and subtract from
#'     the upper bound.
#'
#' @author David Gerard
#'
#' @noRd
lrt_dr_npp_glpknown4 <- function(gl, g1, g2, drbound = 1/6, xi1 = 1/3, xi2 = 1/3, ntry = 5, fudge = 0) {

  ## Same scenario when no 2
  if (g1 != 2  && g2 != 2) {
    return(lrt_dr_pp_glpknown4(gl = gl, g1 = g1, g2 = g2, drbound = drbound, ntry = ntry))
  }

  ## MLE under alternative
  log_qhat1 <- c(em_li(B = gl))
  l1 <- llike_li(B = gl, lpivec = log_qhat1)

  ## MLE under null
  oout <- stats::optim(par = drbound / 2,
                       fn = obj_dr_npp_gl,
                       method = "Brent",
                       lower = fudge,
                       upper = drbound,
                       control = list(fnscale = -1),
                       gl = gl,
                       xi1 = xi1,
                       xi2 = xi2,
                       g1 = g1,
                       g2 = g2)
  l0 <- oout$value
  alpha <- oout$par[[1]]

  ## calculate degrees of freedom
  TOL <- 1e-3
  df_alt <- sum(offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = g1, p2 = g2) > TOL | log_qhat1 > log(TOL)) - 1
  df_null <- !onbound(
    g1 = g1,
    g2 = g2,
    alpha = alpha,
    xi1 = xi1,
    xi2 = xi2,
    dr = TRUE,
    pp = FALSE,
    drbound = drbound,
    TOL = TOL)
  df <- df_alt - df_null
  df <- max(df, 1)

  ## calculate LRT statistic
  llr <- -2 * (l0 - l1)

  ## calculate p-value
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

## pp, no dr ----------------------------------------------------------------

#' Objective for lrt_ndr_pp_glpknown4
#'
#' @inheritParams obj_dr_pp_gl
#' @param par if g1 != 2 or g2 != 2, then should be xi1 (or xi2). If
#'     both g1 == 2 and g2 == 2, then first element is xi1 and second
#'     element is xi2.
#' @param alpha The known rate of double reduction.
#'
#' @author David Gerard
#'
#' @noRd
obj_ndr_pp_gl <- function(par, gl, g1, g2, alpha = 0) {
  par[par < 0] <- 0 ## deal with -1e-18 etc
  par[par > 1] <- 1
  if (g1 != 2 && g2 != 2) {
    stop("obj_ndr_pp_gl: have to have g1 == 2 or g2 == 2")
  } else if (g1 == 2 && g2 != 2) {
    stopifnot(length(par) == 1)
    obj <- like_glpknown_2(
      gl = gl,
      alpha = alpha,
      xi1 = par[[1]],
      xi2 = 1/3,
      g1 = g1,
      g2 = g2,
      log_p = TRUE)
  } else if (g1 != 2 && g2 == 2) {
    stopifnot(length(par) == 1)
    obj <- like_glpknown_2(
      gl = gl,
      alpha = alpha,
      xi1 = 1/3,
      xi2 = par[[1]],
      g1 = g1,
      g2 = g2,
      log_p = TRUE)
  } else if (g1 == 2 && g2 == 2) {
    stopifnot(length(par) == 2)
    obj <- like_glpknown_2(
      gl = gl,
      alpha = alpha,
      xi1 = par[[1]],
      xi2 = par[[2]],
      g1 = g1,
      g2 = g2,
      log_p = TRUE)
  }
  return(obj)
}

#' LRT when there is pp, but no dr
#'
#' @inheritParams lrt_dr_pp_glpknown4
#' @param alpha The known rate of double reduction.
#' @param fudge How much to add to the lower bound and subtract from
#'     the upper bound.
#'
#' @author David Gerard
#'
#' @noRd
lrt_ndr_pp_glpknown4 <- function(gl, g1, g2, alpha = 0, fudge = 0) {
  if (g1 != 2 && g2 != 2) {
    return(lrt_ndr_npp_glpknown4(gl = gl, g1 = g1, g2 = g2, alpha = alpha))
  }

  ## MLE under alternative
  log_qhat1 <- c(em_li(B = gl))
  l1 <- llike_li(B = gl, lpivec = log_qhat1)

  ## max likelihood under null
  if (g1 == 2 && g2 != 2) {
    oout <- stats::optim(
      par = 1/3,
      fn = obj_ndr_pp_gl,
      lower = fudge,
      upper = 1 - fudge,
      method = "Brent",
      control = list(fnscale = -1),
      gl = gl,
      alpha = alpha,
      g1 = g1,
      g2 = g2)
    xi1 <- oout$par[[1]]
    xi2 <- NA_real_
  } else if (g1 != 2 && g2 == 2) {
    oout <- stats::optim(
      par = 1/3,
      fn = obj_ndr_pp_gl,
      lower = fudge,
      upper = 1 - fudge,
      method = "Brent",
      control = list(fnscale = -1),
      gl = gl,
      alpha = alpha,
      g1 = g1,
      g2 = g2)
    xi1 <- NA_real_
    xi2 <- oout$par[[1]]
  } else if (g1 == 2 && g2 == 2) {
    oout <- stats::optim(
      par = c(1/3, 1/3),
      fn = obj_ndr_pp_gl,
      lower = c(fudge, fudge),
      upper = c(1 - fudge, 1 - fudge),
      method = "L-BFGS-B",
      control = list(fnscale = -1),
      gl = gl,
      alpha = alpha,
      g1 = g1,
      g2 = g2)
    xi1 <- oout$par[[1]]
    xi2 <- oout$par[[2]]
  }
  l0 <- oout$value

  ## calculate degrees of freedom
  TOL <- 1e-3
  df_alt <- sum(offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = g1, p2 = g2) > TOL | log_qhat1 > log(TOL)) - 1
  df_null <- !onbound(
    g1 = g1,
    g2 = g2,
    alpha = alpha,
    xi1 = xi1,
    xi2 = xi2,
    dr = FALSE,
    pp = TRUE,
    drbound = 1,
    TOL = TOL)
  df <- df_alt - df_null
  df <- max(df, 1)

  ## calculate LRT statistic
  llr <- -2 * (l0 - l1)

  ## calculate p-value
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
