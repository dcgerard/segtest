#' Iterator over array
#'
#' @param obj An array
#' @param by The dimension to iterate over.
#' @param recycle Should the iterator reset?
#' @param ... not used
#'
#' @author David Gerard
#'
#' @seealso [nextElem.arrayiter()]
#'
#' @examples
#' glist <- multidog_to_g(mout = ufit, type = "all_gl", p1 = "indigocrisp", p2 = "sweetcrisp")
#' g <- iter(glist$g, by = 3)
#' head(nextElem(g))
#' head(nextElem(g))
#' head(nextElem(g))
#'
#' @exportS3Method iterators::iter
iter.array <- function(obj, by = 1, recycle = FALSE, ...) {
  checkFunc <- function(...) TRUE ## hard coding to TRUE
  stopifnot(inherits(obj, "array"))
  stopifnot(by >= 1, by <= length(dim(obj)))
  state <- new.env()
  state$i <- 0L
  state$obj <- obj
  n <- dim(obj)[[by]]
  it <- list(state = state, by = by, length = n, checkFunc = checkFunc, recycle = recycle)
  class(it) <- c("arrayiter", "abstractiter", "iter")
  return(it)
}

#' nextelement in an array over array
#'
#' @param obj An arrayiter object
#' @param ... not used
#'
#' @author David Gerard
#'
#' @seealso [iter.array()]
#'
#' @exportS3Method iterators::nextElem
nextElem.arrayiter <- function(obj, ...) {
  repeat {
    tryCatch({
      obj$state$i <- obj$state$i + 1L
      indlist <- lapply(X = dim(obj$state$obj), FUN = seq_len)
      indlist[[obj$by]] <- obj$state$i
      indlist <- c(list(x = obj$state$obj), indlist)
      return(do.call(what = `[`, args = indlist))
    }, error = function(e) {
            if (any(nzchar(e$message))) {
                if (identical(e$message, "subscript out of bounds")) {
                  if (obj$recycle) {
                    obj$state$i <- 0L
                  }
                  else {
                    stop("StopIteration", call. = FALSE)
                  }
                }
                else {
                  stop(e$message, call. = FALSE)
                }
            }
            else {
                stop("Abort", call. = e)
            }
        })
    }
}


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
#' @param ploidy The ploidy. Note that most methods in this package
#'     (including those in \code{\link{multi_lrt}()}) assume that the ploidy
#'     is 4. But we allow for arbitrary ploidy in this function since it
#'     might be useful in the future.
#'
#' @return A list with the following elements
#' \describe{
#'    \item{g}{Either a matrix of counts, where the columns index the genotype
#'             and the rows index the loci (\code{type = "all_g"} or
#'             \code{type = "off_g"}). Or an array of genotype (natural) log-likelihoods
#'             where the rows index the loci, the columns index the
#'             individuals, and the slices index the genotypes
#'             (\code{type = "all_gl"} or \code{type = "off_gl"}).}
#'    \item{p1}{Either a vector of known parental genotypes
#'              (\code{type = "off_gl"}, \code{type = "all_g"} or \code{type = "off_g"}).
#'              Or a matrix of genotype (natural) log-likelihoods where the
#'              rows index the loci and the columns index the genotypes
#'              (\code{type = "all_gl"}).}
#'    \item{p2}{Either a vector of known parental genotypes
#'              (\code{type = "off_gl"}, \code{type = "all_g"} or \code{type = "off_g"}).
#'              Or a matrix of genotype (natural) log-likelihoods where the
#'              rows index the loci and the columns index the genotypes
#'              (\code{type = "all_gl"}).}
#' }
#'
#' @author David Gerard
#'
#' @examples
#' multidog_to_g(mout = ufit, type = "all_g", p1 = "indigocrisp", p2 = "sweetcrisp")
#' multidog_to_g(mout = ufit, type = "all_gl", p1 = "indigocrisp", p2 = "sweetcrisp")
#'
#'
#' @export
multidog_to_g <- function(
    mout,
    type = c("off_gl", "all_gl", "all_g", "off_g"),
    p1 = NULL,
    p2 = NULL,
    ploidy = 4) {
  type <- match.arg(type)

  if (type == "off_gl") {
    stopifnot(!is.null(mout$snpdf$p1geno), !is.null(mout$snpdf$p2geno))
    p1_geno <- mout$snpdf$p1geno
    p2_geno <- mout$snpdf$p2geno
    g <- updog::format_multidog(ufit, varname = paste0("logL_", 0:ploidy))
  } else if (type == "all_gl") {
    stopifnot(!is.null(p1), !is.null(p2))
    g <- updog::format_multidog(ufit, varname = paste0("logL_", 0:ploidy))
    p1_geno <- g[, p1, ]
    p2_geno <- g[, p2, ]
    g <- g[, !(dimnames(g)[[2]] %in% c(p1, p2)), ]
  } else if (type == "off_g") {
    stopifnot(!is.null(mout$snpdf$p1geno), !is.null(mout$snpdf$p2geno))
    p1_geno <- mout$snpdf$p1geno
    p2_geno <- mout$snpdf$p2geno
    gmat <- updog::format_multidog(ufit, varname = "geno")
    g <- t(apply(X = gmat, MARGIN = 1, FUN = gvec_to_gcount, ploidy = ploidy))
    colnames(g) <- 0:ploidy
  } else if (type == "all_g") {
    stopifnot(!is.null(p1), !is.null(p2))
    gmat <- updog::format_multidog(ufit, varname = "geno")
    p1_geno <- gmat[, p1, drop = TRUE]
    p2_geno <- gmat[, p2, drop = TRUE]
    gmat <- gmat[, !(colnames(gmat) %in% c(p1, p2)), drop = FALSE]
    g <- t(apply(X = gmat, MARGIN = 1, FUN = gvec_to_gcount, ploidy = ploidy))
    colnames(g) <- 0:ploidy
  }

  return(list(g = g, p1 = p1_geno, p2 = p2_geno))
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
#' glist <- multidog_to_g(mout = ufit, type = "all_g", p1 = "indigocrisp", p2 = "sweetcrisp")
#' p1 <- glist$p1
#' p2 <- glist$p2
#' g <- glist$g
#'
#' ## Using genotype likelihoods (typically good idea)
#' glist <- multidog_to_g(mout = ufit, type = "all_gl", p1 = "indigocrisp", p2 = "sweetcrisp")
#' p1 <- glist$p1
#' p2 <- glist$p2
#' g <- glist$g
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

  if (type == "g") {
    ret <- foreach::foreach(
      g_now = iter(g, by = "row"),
      p1_now = iter(p1),
      p2_now = iter(p2),
      .combine = rbind
    ) %dorng% {
        lout <- lrt_men_g4(
          x = g_now,
          g1 = p1_now,
          g2 = p2_now,
          drbound = drbound,
          pp = pp,
          dr = dr,
          alpha = alpha,
          xi1 = xi1,
          xi2 = xi2)
        unlist(lout)
    }
  } else if (type == "gl") {
    ret <- foreach::foreach(
      g_now = iter(g, by = 1),
      p1_now = iter(p1, by = "row"),
      p2_now = iter(p2, by = "row"),
      .combine = rbind
    ) %dorng% {
        lout <- lrt_men_gl4(
          gl = g_now,
          g1 = p1_now,
          g2 = p2_now,
          drbound = drbound,
          pp = pp,
          dr = dr,
          alpha = alpha,
          xi1 = xi1,
          xi2 = xi2)
        unlist(lout)
    }
  }

  attr(ret, "rng") <- NULL
  attr(ret, "doRNG_version") <- NULL

  return(ret)
}
