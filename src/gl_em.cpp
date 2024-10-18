#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//' Log-sum-exponential trick.
//'
//' @param x A vector to log-sum-exp.
//'
//' @return The log of the sum of the exponential
//'     of the elements in \code{x}.
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
double log_sum_exp(const arma::vec &x) {
  double max_x = arma::max(x);
  double lse; // the log-sum-exp
  // if all -Inf, need to treat this special to avoid -Inf + Inf.
  if (max_x == -arma::datum::inf) {
    lse = -arma::datum::inf;
  } else {
    lse = max_x + std::log(arma::sum(arma::exp(x - max_x)));
  }
  return lse;
}

//' Log-sum-exponential trick using just two doubles.
//'
//' @param x A double.
//' @param y Another double.
//'
//' @return The log of the sum of the exponential of x and y.
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
double log_sum_exp_2(double x, double y) {
  double z = std::max(x, y);
  double finalval;
  if (z == -arma::datum::inf) {
    finalval = -arma::datum::inf;
  } else {
    finalval = std::log(std::exp(x - z) + std::exp(y - z)) + z;
  }
  return finalval;
}

//' Objective function for \code{\link{em_li}()}
//'
//' @param B The log-likelihood matrix. Rows are individuals columns are
//'     genotypes.
//' @param lpivec The log prior vector.
//'
//' @author David Gerard
//'
//' @return The log-likelihood of a vector of genotype frequencies when
//'     using genotype likelihoods. This is from Li (2011).
//'
//' @references
//' \itemize{
//'   \item{Li, H. (2011). A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. \emph{Bioinformatics}, 27(21), 2987-2993. \doi{10.1093/bioinformatics/btr509}}
//' }
//'
//' @examples
//' # Simulate some data
//' set.seed(1)
//' gl <- simgl(nvec = c(3, 2, 4, 1, 2))
//' # Log-likelihood at given log-priors
//' prob <- c(0.1, 0.2, 0.4, 0.2, 0.1)
//' lprob <- log(prob)
//' llike_li(B = gl, lpivec = lprob)
//'
//' @export
// [[Rcpp::export]]
double llike_li(const arma::mat &B, const arma::vec &lpivec) {

  double ll = 0.0;
  int K = B.n_cols - 1;
  int n = B.n_rows;

  if (B.n_cols != lpivec.n_elem) {
    Rcpp::stop("Number of columns in B should equal length of lpivec");
  }

  for (int i = 0; i < n; i++) {
    double ival = R_NegInf;
    for (int j = 0; j <= K; j++) {
      ival = log_sum_exp_2(ival, B(i, j) + lpivec(j));
    }
    ll += ival;
  }
  return ll;
}

//' EM algorithm from Li (2011)
//'
//' EM algorithm to estimate prior genotype probabilities from genotype
//' likelihoods.
//'
//' @param B Matrix of genotype log-likelihoods. The rows index the individuals
//'     and the columns index the genotypes.
//' @param itermax The maximum number of iterations.
//' @param eps The stopping criteria.
//'
//' @return A vector of log prior probabilities for each genotype.
//'
//' @author David Gerard
//'
//' @references
//' \itemize{
//'   \item{Li, H. (2011). A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. \emph{Bioinformatics}, 27(21), 2987-2993. \doi{10.1093/bioinformatics/btr509}}
//' }
//'
//' @examples
//' # Simulate some data
//' set.seed(1)
//' gl <- simgl(nvec = c(3, 2, 4, 1, 2))
//' # Run em
//' lprob <- em_li(B = gl)
//' # Exponentiate to get probabilities
//' prob <- exp(c(lprob))
//' prob
//'
//' @export
// [[Rcpp::export]]
arma::vec em_li(const arma::mat &B, int itermax = 100, double eps = 1e-5) {
  int K = B.n_cols - 1;
  int n = B.n_rows;
  double valinit = -std::log((double)K + 1);

  arma::vec lpivec(K + 1, arma::fill::value(valinit));

  double llold;
  double llnew = llike_li(B, lpivec);
  double err = R_PosInf;
  arma::vec lw(K + 1);

  int i = 0;
  while((i < itermax) & (err > eps)) {
    llold = llnew;
    lw.fill(-arma::datum::inf);

    for (int j = 0; j < n; j++) {
      arma::vec api(K + 1);
      for (int k = 0; k <= K; k++) {
        api(k) = B(j, k) + lpivec(k);
      }
      api = api - log_sum_exp(api);
      for (int k = 0; k <= K; k++) {
        lw(k) = log_sum_exp_2(lw(k), api(k));
      }
    }

    lpivec = lw - log_sum_exp(lw);
    llnew = llike_li(B, lpivec);
    if (llnew < llold) {
      Rcpp::stop("log-likelihood not increasing");
    }

    err = std::abs(llnew - llold);
    i++;
  }

  return lpivec;
}
