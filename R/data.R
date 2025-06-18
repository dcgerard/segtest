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

#' Disomic and polysomic segregation patterns
#'
#' Gamete frequencies for all possible disomic and polysomic segregation
#' patterns for even ploidies 2 through 20. If you need higher ploidy levels,
#' let me know and I'll update it (it's very easy).
#'
#' @format A data frame with the following columns
#' \describe{
#'   \item{\code{ploidy}}{The ploidy of the parent.}
#'   \item{\code{g}}{The genotype of the parent.}
#'   \item{\code{m}}{The pairing configuration given disomic inheritance. See Gerard et al (2018).}
#'   \item{\code{p}}{The gamete frequencies. Element \code{p[[i]]} is the probability a gamete will have dosage \code{i-1}.}
#'   \item{\code{mode}}{Whether the inheritance pattern that leads to these gamete frequencies is \code{"disomic"}, \code{"polysomic"}, or \code{"both"}.}
#' }
#'
#' @references
#' \itemize{
#'   \item{Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018). Genotyping polyploids from messy sequencing data. Genetics, 210(3), 789-807.}
#' }
#'
#' @author David Gerard
"seg"
