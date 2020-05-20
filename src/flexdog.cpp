#include "mupdog.h"

//' E-step in \code{\link{flexdog}}.
//'
//' @inheritParams flexdog_full
//' @param probk_vec The vector of current prior probabilities of each genotype.
//'
//'
//' @return A matrix of numerics. The rows index the individuals and the
//'     columns index the genotype. These weights are used in the EM algorithm
//'     (and is indeed the E-step) in \code{\link{flexdog_full}}.
//'
//' @author David Gerard
//'
//' @keywords internal
//' @noRd
//'
//' @seealso \code{\link{flexdog}} for the full EM algorithm.
//'
// [[Rcpp::export]]
NumericMatrix get_wik_mat(NumericVector probk_vec,
                          NumericVector refvec,
                          NumericVector sizevec,
                          int ploidy,
                          double seq,
                          double bias,
                          double od) {
  // Check input -----------------------------------------------------------
  int nind = refvec.length();
  if (nind != sizevec.length()) {
    Rcpp::stop("get_wik_mat: sizevec and refvec must have the same length.");
  }
  if (probk_vec.length() != ploidy + 1) {
    Rcpp::stop("get_wik_mat: probk_vec must have length ploidy + 1.");
  }

  // Calculate the posterior probability of each genotype -------------------
  NumericMatrix wik_mat(nind, ploidy + 1);
  NumericVector lprobk_vec = Rcpp::log(probk_vec);
  NumericVector xi(ploidy + 1);
  for (int k = 0; k <= ploidy; k++) {
    xi(k) = xi_double((double)k / (double)ploidy, seq, bias);
  }

  double sumi; // denominator to get wik for each i.
  NumericVector wvec(ploidy + 1);
  for (int i = 0; i < nind; i++) {
    for (int k = 0; k <= ploidy; k++) {
      wvec(k) = lprobk_vec(k) + dbetabinom_double(refvec(i), sizevec(i), xi(k), od, true);
    }
    sumi = log_sum_exp(wvec);
    wvec = Rcpp::exp(wvec - sumi);
    wik_mat(i, _) = wvec;
  }

  return wik_mat;
}

//' Returns genotype log-likelihood matrix
//'
//' @inheritParams get_wik_mat
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
NumericMatrix get_genotype_likelihoods(NumericVector refvec,
                                       NumericVector sizevec,
                                       int ploidy,
                                       double seq,
                                       double bias,
                                       double od) {
  // Check input -----------------------------------------------------------
  int nind = refvec.length();
  if (nind != sizevec.length()) {
    Rcpp::stop("get_wik_mat: sizevec and refvec must have the same length.");
  }

  // Calculate the posterior probability of each genotype -------------------
  NumericMatrix wik_mat(nind, ploidy + 1);
  NumericVector xi(ploidy + 1);
  for (int k = 0; k <= ploidy; k++) {
    xi(k) = xi_double((double)k / (double)ploidy, seq, bias);
  }

  for (int i = 0; i < nind; i++) {
    for (int k = 0; k <= ploidy; k++) {
      wik_mat(i, k) = dbetabinom_double(refvec(i), sizevec(i), xi(k), od, true);
    }
  }

  return wik_mat;
}

//' Log-likelihood that \code{\link{flexdog}} maximizes.
//'
//' @inheritParams flexdog_full
//' @param probk_vec The kth element is the prior probability of genotype k (when starting to count from 0).
//'
//' @author David Gerard
//'
//' @keywords internal
//' @noRd
//'
//' @return The objective (marginal log-likelihood) used in
//'     \code{\link{flexdog_full}}.
//'
// [[Rcpp::export]]
double flexdog_obj(NumericVector probk_vec,
                   NumericVector refvec,
                   NumericVector sizevec,
                   int ploidy,
                   double seq,
                   double bias,
                   double od,
                   double mean_bias,
                   double var_bias,
                   double mean_seq,
                   double var_seq,
                   double mean_od,
                   double var_od) {
  // Check input -----------------------------------------------------------
  int nind = refvec.length();
  if (nind != sizevec.length()) {
    Rcpp::stop("get_wik_mat: sizevec and refvec must have the same length.");
  }
  if (probk_vec.length() != ploidy + 1) {
    Rcpp::stop("get_wik_mat: probk_vec must have length ploidy + 1.");
  }

  // Calculate the posterior probability of each genotype -------------------
  NumericVector lprobk_vec = Rcpp::log(probk_vec);
  NumericVector xi(ploidy + 1);
  for (int k = 0; k <= ploidy; k++) {
    xi(k) = xi_double((double)k / (double)ploidy, seq, bias);
  }

  // Calculate likelihood ---------------------------------------------------
  double obj = 0.0;
  NumericVector wvec(ploidy + 1);
  for (int i = 0; i < nind; i++) {
    for (int k = 0; k <= ploidy; k++) {
      wvec(k) = lprobk_vec(k) + dbetabinom_double(refvec(i), sizevec(i), xi(k), od, true);
    }
    obj = obj + log_sum_exp(wvec);
  }

  // Penalties --------------------------------------------------------------
  obj = obj + pen_bias(bias, mean_bias, var_bias);
  obj = obj + pen_seq_error(seq, mean_seq, var_seq);
  obj = obj + pen_seq_error(od, mean_od, var_od);
  return obj;
}

//' Objective for mixture of known dist and uniform dist.
//'
//' @param alpha The mixing weight.
//' @param pvec The known distribution (e.g. from assuming an
//'     F1 population).
//' @param weight_vec A vector of weights.
//'
//' @return The objective when updating \code{pivec} when \code{model = "f1"}
//'     or \code{model = "s1"} in \code{\link{flexdog_full}}.
//'
//' @author David Gerard
//'
//' @keywords internal
//' @noRd
//'
// [[Rcpp::export]]
double f1_obj(double alpha,
              arma::vec pvec,
              arma::vec weight_vec) {

  // check input ----------------------------------
  int ploidy = pvec.n_elem - 1;
  if (weight_vec.n_elem != (unsigned)ploidy + 1) {
    Rcpp::stop("f1_obj: pvec and weight_vec should be the same length.");
  }
  if ((alpha < 0) | (alpha > 1)) {
    Rcpp::stop("f1_obj: alpha should be between 0 and 1.");
  }

  // get obj --------------------------------------
  double obj = 0.0;
  double newp;
  for (int k = 0; k <= ploidy; k++) {
    newp = (1.0 - alpha) * pvec(k) + alpha / ((double)ploidy + 1.0);
    obj = obj + weight_vec(k) * std::log(newp);
  }
  return obj;
}
