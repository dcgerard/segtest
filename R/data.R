#' Genotype data from Cappai et al. (2020)
#'
#' A subset of data from Cappai et al. (2020), fit using
#' \code{\link[updog]{multidog}()}. This just contains a random set of 10 loci.
#'
#' @format An object of type \code{multidog} output from \code{\link[updog]{multidog}()}.
#' \describe{
#'   \item{\code{ufit}}{Uses the \code{model = "norm"} option.}
#'   \item{\code{ufit2}}{Uses the \code{model = "f1pp"} option.}
#'   \item{\code{ufit3}}{Uses the \code{model = "f1"} option.}
#' }
#'
#' @source \doi{10.5281/zenodo.13715703}
#'
#' @references
#' \itemize{
#'   \item{Cappai, F., Amadeu, R. R., Benevenuto, J., Cullen, R., Garcia, A., Grossman, A., Ferrão, L., & Munoz, P. (2020). High-resolution linkage map and QTL analyses of fruit firmness in autotetraploid blueberry. \emph{Frontiers in plant science}, 11, 562171. \doi{10.3389/fpls.2020.562171}.}
#' }
"ufit"

#' @rdname ufit
"ufit2"

#' @rdname ufit
"ufit3"
