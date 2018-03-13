#include "mupdog.h"

// functions for fitting flexdog


//' Obtain the genotype distribution given the distribution of discrete uniforms.
//'
//' @inheritParams flexdog
//' @param pivec The mixing probability of the i'th discrete uniform distribution.
//'
//' @author David Gerard
//'
//' @seealso \code{\link{flexdog}} where this is used.
//'
// [[Rcpp::export]]
NumericVector get_probk_vec(NumericVector pivec, std::string model, double mode) {
  int K = pivec.length() - 1;
  NumericVector probk_vec(K + 1);
  if ((model == "flex") | (model == "hw")) {
    probk_vec = pivec;
  } else if (model == "ash") {
    double denom; // what you divide the pi's by.
    for (int i = 0; i <= K; i++) { // iterate through pivec
      if (abs((double)i - mode) < TOL) {
        probk_vec(i) += pivec(i);
      } else if (i < mode) {
        denom = (double)(std::floor(mode) - i + 1);
        for (int j = i; j <= mode; j++) { //iterate through probk_vec
          probk_vec(j) += pivec(i) / denom;
        }
      } else if (i > mode) { //iterate through probk_vec
        denom = (double)(i - std::ceil(mode) + 1);
        for (int j = i; j >= mode; j--) {
          probk_vec(j) += pivec(i) / denom;
        }
      } else {
        Rcpp::stop("get_probk_vec: How did you get here??");
      }
    }
  } else {
    Rcpp::stop("get_probk_vec: model must be one of 'flex', 'ash', or 'hw'");
  }
  return probk_vec;
}

//' Compute inner weights for updating the mixing proportions when using ash model.
//'
//' The (i,k)th element is \eqn{1(k \in F(a, i)) / |F(a,i)|}.
//'
//' @inheritParams flexdog
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
NumericMatrix get_inner_weights(int ploidy, double mode) {
  NumericMatrix inner_weights(ploidy + 1, ploidy + 1);
  double denom;
  for (int i = 0; i <= ploidy; i++) {
    if (std::abs((double)i - mode) < TOL) {
      inner_weights(i, i) = 1.0;
    } else if (i < mode) {
      denom = (double)(std::floor(mode) - i + 1);
      for (int j = i; j <= mode; j++) {
        inner_weights(i, j) = 1.0 / denom;
      }
    } else if (i > mode) {
      denom = (double)(i - std::ceil(mode) + 1);
      for (int j = i; j >= mode; j--) {
        inner_weights(i, j) = 1.0 / denom;
      }
    } else {
      Rcpp::stop("get_inner_weights: How did you get here??");
    }
  }
  return(inner_weights);
}


//' E-step in \code{\link{flexdog}}.
//'
//' @inheritParams flexdog
//' @param probk_vec The vector of current prior probabilities of each genotype.
//'
//'
//' @author David Gerard
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

//' Log-likelihood that \code{\link{flexdog}} maximizes.
//'
//' @inheritParams flexdog
//' @param probk_vec The kth element is the prior probability of genotype k (when starting to count from 0).
//'
//' @author David Gerard
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
                   double var_seq) {
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
  return obj;
}



















