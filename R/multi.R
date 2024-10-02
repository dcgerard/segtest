#' Converts multidog output to a format usable for multi_lrt()
#'
#' @param mout The output of \code{\link[updog]{multidog}()}.
#' @param type
#' \describe{
#'   \item{\code{"off_gl"}}{Genotype likelihoods of offspring but not parents. This is the typical choice if you used the "f1" or "f1pp" options when genotyping.}
#'   \item{\code{"all_gl"}}{Genotype likelihoods of offspring and parents. This is only done if you did \emph{not} use the "f1" or "f1pp" options when genotyping. If this is the case, then you need to specify which individuals are the parents.}
#'   \item{\code{"off_g"}}{Genotypes, assuming that they are known. You used the "f1" or "f1pp" option when genotyping.}
#'   \item{\code{"all_g"}}{Genotypes, assuming that they are known. You did \emph{not} use the "f1" or "f1pp" option when genotyping. If this is the case, then you need to specify which individuals are the parents}.
#' }
#' @param p1 The first parent name if using \code{type = "all_gl"} or \code{type = "all_g"}.
#' @param p1 The second parent name if using \code{type = "all_gl"} or \code{type = "all_g"}.
#'
#' @author David Gerard
#'
#' @export
multidog_to_g <- function(
    mout,
    type = c("off_gl", "all_gl", "all_g", "off_g"),
    p1 = NULL,
    p2 = NULL) {
  type <- match.arg(type)

  if (type == "off_gl") {

  } else if (type == "all_gl") {

  } else if (type == "off_g") {

  } else if (type == "all_g") {

  }

}

#' Parallelized likelihood ratio test for segregation distortion.
#'
#' Uses the \code{future} package to implement parallelization support for
#' the likelihood ratio tests for segregation distortion. Right now, this is
#' only supported for tetraploids (allo, auto, or segemental).
#'
#' @param g One of two inputs
#'   \itemize{
#'     \item{A matrix of genotype counts. The rows index the locis and the columns index the genotypes.}
#'     \item{An array of genotype log-likelihoods. The rows index the loci, the columns index the individuals, and the slices index the genotypes. Log-likelihoods are base e (natural log).}
#'   }
#' @param p1 One of two inputs
#'   \itemize{
#'     \item{A vector of parent 1's genotypes.}
#'     \item{A matrix of parent 1's genotype log-likelihoods. The rows index the loci and the columns index the genotypes. Logs are in base e (natural log).}
#'   }
#' @param p2 One of two inputs
#'   \itemize{
#'     \item{A vector of parent 1's genotypes.}
#'     \item{A matrix of parent 1's genotype log-likelihoods. The rows index the loci and the columns index the genotypes. Logs are in base e (natural log).}
#'   }
#' @inheritParams lrt_men_gl4
#'
#' @author David Gerard
#'
#' @examples
#' ## Assuming genotypes are known (typically bad idea)
#' gmat <- updog::format_multidog(ufit, varname = "geno")
#' p1 <- gmat[, "indigocrisp"]
#' p2 <- gmat[, "sweetcrisp"]
#' gmat <- gmat[, !(colnames(g) %in% c("indigocrisp", "sweetcrisp"))]
#'
#' ## Using genotype likelihoods (typically good idea)
#' g <- updog::format_multidog(ufit, varname = paste0("logL_", 0:4))
#' p1 <- g[, "indigocrisp", ]
#' p2 <- g[, "sweetcrisp", ]
#' g <- g[, !(dimnames(g)[[2]] %in% c("indigocrisp", "sweetcrisp")), ]
#'
#'
#' @export
multi_lrt <- function(g,
                      p1,
                      p2,
                      drbound = 1/6,
                      pp = TRUE,
                      dr = TRUE,
                      alpha = 0,
                      xi1 = 1/3,
                      xi2 = 1/3) {
  ## Check input --------------------------------------------------------------
  if (length(dim(g)) == 3) {
    type <- "gl"
    nloc <- dim(g)[[1]]
    nind <- dim(g)[[2]]
    stopifnot(dim(g)[[3]] == 5)
  } else if (length(dim(g)) == 2) {
    type <- "g"
    nloc <- nrow(g)
    stopifnot(ncol(g) == 5)
  }

  ## Register doFuture  -------------------------------------------------------
  oldDoPar <- doFuture::registerDoFuture()
  on.exit(with(oldDoPar, foreach::setDoPar(fun=fun, data=data, info=info)), add = TRUE)

  ret <- foreach::foreach(
    i = seq_len(nloc)
  ) %dopar% {
    if (type == g) {
      lout <- lrt_men_g4(
        x = g[i, ],
        g1 = p1[[i]],
        g2 = p2[[i]],
        drbound = drbound,
        pp = pp,
        dr = dr,
        alpha = alpha,
        xi1 = xi1,
        xi2 = xi2)
    }
  }

}
