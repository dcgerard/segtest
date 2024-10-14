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
#' g <- iterators::iter(glist$g, by = 3)
#' head(iterators::nextElem(g))
#' head(iterators::nextElem(g))
#' head(iterators::nextElem(g))
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

#' Next element in an array
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
#' @param p2 The second parent name if using \code{type = "all_gl"} or \code{type = "all_g"}.
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
#' multidog_to_g(mout = ufit2, type = "off_g")
#' multidog_to_g(mout = ufit2, type = "off_gl")
#' multidog_to_g(mout = ufit3, type = "off_g")
#' multidog_to_g(mout = ufit3, type = "off_gl")
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
    if (!is.null(mout$snpdf$p1geno) && !is.null(mout$snpdf$p2geno)) {
      p1_geno <- mout$snpdf$p1geno
      p2_geno <- mout$snpdf$p2geno
    } else if (!is.null(mout$snpdf$ell1) && !is.null(mout$snpdf$ell2)) {
      p1_geno <- mout$snpdf$ell1
      p2_geno <- mout$snpdf$ell2
    } else {
      stop("mout was not fit using either the 'f1' or the 'f1pp' models")
    }
    g <- updog::format_multidog(mout, varname = paste0("logL_", 0:ploidy))
  } else if (type == "all_gl") {
    stopifnot(!is.null(p1), !is.null(p2))
    g <- updog::format_multidog(mout, varname = paste0("logL_", 0:ploidy))
    p1_geno <- g[, p1, ]
    p2_geno <- g[, p2, ]
    g <- g[, !(dimnames(g)[[2]] %in% c(p1, p2)), ]
  } else if (type == "off_g") {
    if (!is.null(mout$snpdf$p1geno) && !is.null(mout$snpdf$p2geno)) {
      p1_geno <- mout$snpdf$p1geno
      p2_geno <- mout$snpdf$p2geno
    } else if (!is.null(mout$snpdf$ell1) && !is.null(mout$snpdf$ell2)) {
      p1_geno <- mout$snpdf$ell1
      p2_geno <- mout$snpdf$ell2
    } else {
      stop("mout was not fit using either the 'f1' or the 'f1pp' models")
    }
    gmat <- updog::format_multidog(mout, varname = "geno")
    g <- t(apply(X = gmat, MARGIN = 1, FUN = gvec_to_gcount, ploidy = ploidy))
    colnames(g) <- 0:ploidy
  } else if (type == "all_g") {
    stopifnot(!is.null(p1), !is.null(p2))
    gmat <- updog::format_multidog(mout, varname = "geno")
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
#' only supported for tetraploids (allo, auto, or segmental). This function
#' is only somewhat tested (the single-locus LRT functions in the "See Also"
#' section are very well tested). So please send any bugs you notice to
#' \url{https://github.com/dcgerard/segtest/issues}.
#'
#' @section Parallel Computation:
#'
#' The \code{multi_lrt()} function supports parallel computing. It does
#' so through the \href{https://cran.r-project.org/package=future}{future}
#' package.
#'
#' You first specify the evaluation plan with \code{\link[future]{plan}()}
#' from the \code{future} package. On a local machine, this is typically
#' just \code{future::plan(future::multisession, workers = nc)} where
#' \code{nc} is the number of workers you want. You can find the maximum
#' number of possible workers with \code{\link[future]{availableCores}()}.
#' You then run \code{multi_lrt()}, then shut down the workers with
#' \code{future::plan(future::sequential)}.
#'
#' @param g One of two inputs
#'   \itemize{
#'     \item{A matrix of genotype counts. The rows index the loci and the columns index the genotypes.}
#'     \item{An array of genotype log-likelihoods. The rows index the loci, the columns index the individuals, and the slices index the genotypes. Log-likelihoods are base e (natural log).}
#'   }
#' @param p1 One of three inputs
#'   \itemize{
#'     \item{A vector of parent 1's genotypes.}
#'     \item{A matrix of parent 1's genotype log-likelihoods. The rows index the loci and the columns index the genotypes. Logs are in base e (natural log).}
#'     \item{\code{NULL} (only supported when using genotype likelihoods for the offspring)}
#'   }
#' @param p2 One of three inputs
#'   \itemize{
#'     \item{A vector of parent 1's genotypes.}
#'     \item{A matrix of parent 1's genotype log-likelihoods. The rows index the loci and the columns index the genotypes. Logs are in base e (natural log).}
#'     \item{\code{NULL} (only supported when using genotype likelihoods for the offspring)}
#'   }
#' @param nullprop Should we return the null proportions (\code{TRUE}) or not (\code{FALSE})?
#' @inheritParams lrt_men_gl4
#'
#' @author David Gerard
#'
#' @examples
#' ## Assuming genotypes are known (typically a bad idea)
#' glist <- multidog_to_g(mout = ufit, type = "all_g", p1 = "indigocrisp", p2 = "sweetcrisp")
#' p1_1 <- glist$p1
#' p2_1 <- glist$p2
#' g_1 <- glist$g
#' multi_lrt(g = g_1, p1 = p1_1, p2 = p2_1)
#'
#' ## Using genotype likelihoods (typically a good idea)
#' glist <- multidog_to_g(mout = ufit, type = "all_gl", p1 = "indigocrisp", p2 = "sweetcrisp")
#' p1_2 <- glist$p1
#' p2_2 <- glist$p2
#' g_2 <- glist$g
#' multi_lrt(g = g_2, p1 = p1_2, p2 = p2_2)
#'
#' ## Offspring genotype likelihoods and parent genotypes known
#' multi_lrt(g = g_2, p1 = p1_1, p2 = p2_1)
#'
#' ## Offspring genotype likelihoods and no information on parent genotypes
#' multi_lrt(g = g_2, p1 = NULL, p2 = NULL)
#'
#' \dontrun{
#' ## Parallel computing is supported through the future package
#' future::plan(future::multisession, workers = 2)
#' multi_lrt(g = g_2, p1 = p1_2, p2 = p2_2)
#' future::plan(future::sequential)
#' }
#'
#' @seealso
#' - [lrt_men_g4()] Single locus LRT for segregation distortion when genotypes are known.
#' - [lrt_men_gl4()] Single locus LRT for segregation distortion when using genotype likelihoods.
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
                      xi2 = 1/3,
                      nullprop = FALSE) {

  if (is.null(p1)) {
    p1 <- rep(NA_real_, length.out = dim(g)[[1]])
    names(p1) <- dimnames(g)[[1]]
  }
  if (is.null(p2)) {
    p2 <- rep(NA_real_, length.out = dim(g)[[1]])
    names(p2) <- dimnames(g)[[1]]
  }

  ## Check input --------------------------------------------------------------
  if (length(dim(g)) == 3) {
    if (inherits(p1, "matrix") && inherits(p2, "matrix")) {
      type <- "glp"
      stopifnot(dimnames(g)[[1]] == rownames(p1))
      stopifnot(dimnames(g)[[1]] == rownames(p2))
    } else {
      type <- "glo"
      stopifnot(dimnames(g)[[1]] == names(p1))
      stopifnot(dimnames(g)[[1]] == names(p2))
    }
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

  g_now <- NULL
  p1_now <- NULL
  p2_now <- NULL
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
        as.data.frame(lout)
    }
    if (is.null(ret$p1)) {
      ret <- cbind(ret, p1 = p1)
    }
    if (is.null(ret$p2)) {
      ret <- cbind(ret, p2 = p2)
    }
  } else if (type == "glp") {
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
        as.data.frame(lout)
    }
  } else if (type == "glo") {
    ret <- foreach::foreach(
      g_now = iter(g, by = 1),
      p1_now = iter(p1),
      p2_now = iter(p2),
      .combine = rbind
    ) %dorng% {
        lout <- lrt_men_gl4(
          gl = g_now,
          g1 = if(is.na(p1_now)) NULL else p1_now,
          g2 = if(is.na(p2_now)) NULL else p2_now,
          drbound = drbound,
          pp = pp,
          dr = dr,
          alpha = alpha,
          xi1 = xi1,
          xi2 = xi2)
        as.data.frame(lout)
    }
    if (is.null(ret$p1)) {
      ret <- cbind(ret, p1 = p1)
    }
    if (is.null(ret$p2)) {
      ret <- cbind(ret, p2 = p2)
    }
  }

  ret$snp <- dimnames(g)[[1]]
  rownames(ret) <- NULL

  if (nullprop) {
    ## Add null proportions
    ret$e_pr0 <- NA_real_
    ret$e_pr1 <- NA_real_
    ret$e_pr2 <- NA_real_
    ret$e_pr3 <- NA_real_
    ret$e_pr4 <- NA_real_
    for (i in seq_len(nrow(ret))) {
      tryCatch({
        eval <- offspring_gf_2(alpha = ret$alpha[[i]], xi1 = ret$xi1[[i]], xi2 = ret$xi2[[i]], p1 = ret$p1[[i]], p2 = ret$p2[[i]])
        ret$e_pr0[[i]] <- eval[[1]]
        ret$e_pr1[[i]] <- eval[[2]]
        ret$e_pr2[[i]] <- eval[[3]]
        ret$e_pr3[[i]] <- eval[[4]]
        ret$e_pr4[[i]] <- eval[[5]]
      }, error = function(e){})
    }
  }

  attr(ret, "rng") <- NULL
  attr(ret, "doRNG_version") <- NULL

  return(ret)
}
