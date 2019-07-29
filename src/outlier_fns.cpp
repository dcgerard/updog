#include "mupdog.h"

NumericMatrix get_wik_mat(NumericVector probk_vec,
                          NumericVector refvec,
                          NumericVector sizevec,
                          int ploidy,
                          double seq,
                          double bias,
                          double od);

//' The outlier distribution we use. Right now it is just a
//' beta binomial with mean 1/2 and od 1/3 (so underlying beta
//' is just a uniform from 0 to 1).
//'
//' @param x The number of reference counts.
//' @param n The total number of read-counts.
//' @param logp Return the log density \code{TRUE} or not \code{FALSE}?
//'
//' @return A double. The outlier density value.
//'
//'
//' @author David Gerard
//'
//'
// [[Rcpp::export]]
double doutdist(int x, int n, bool logp) {
  double dval = dbetabinom_double(x, n, 0.5, 1.0 / 3.0, logp);
  return dval;
}

//' E-step in \code{\link{flexdog}} where we now allow an outlier distribution.
//'
//' @inheritParams flexdog_full
//' @param probk_vec The vector of current prior probabilities of each genotype.
//' @param out_prop The probability of being an outlier.
//'
//' @author David Gerard
//'
//' @seealso \code{\link{flexdog}} for the full EM algorithm.
//'     \code{\link{get_wik_mat}} for the equivalent function
//'     without outliers. \code{\link{doutdist}} for the outlier
//'     density function.
//'
//' @return Same as \code{\link{get_wik_mat}} but the last column is
//'     for the outlier class.
//'
// [[Rcpp::export]]
NumericMatrix get_wik_mat_out(NumericVector probk_vec,
                              double out_prop,
                              NumericVector refvec,
                              NumericVector sizevec,
                              int ploidy,
                              double seq,
                              double bias,
                              double od) {
  // Check input -----------------------------------------------------------
  int nind = refvec.length();
  if (nind != sizevec.length()) {
    Rcpp::stop("get_wik_mat_out: sizevec and refvec must have the same length.");
  }
  if (probk_vec.length() != ploidy + 1) {
    Rcpp::stop("get_wik_mat_out: probk_vec must have length ploidy + 1.");
  }


  NumericMatrix wik_mat(nind, ploidy + 2);
  if (out_prop < TOL) { // outliers have zero probability
    NumericMatrix wik_sub(nind, ploidy + 1);
    wik_sub = get_wik_mat(probk_vec,
                          refvec,
                          sizevec,
                          ploidy,
                          seq,
                          bias,
                          od);

    for (int i = 0; i < nind; i++) {
      for (int k = 0; k <= ploidy; k++) {
        wik_mat(i, k) = wik_sub(i, k);
      }
    }
  } else if (1.0 - out_prop < TOL) {
    for (int i = 0; i < nind; i++) { // Last column is 1, rest are 0 since all are outliers.
      wik_mat(i, ploidy + 1) = 1.0;
    }
  } else {
    // Calculate the posterior probability of each genotype -------------------
    NumericVector lprobk_vec = Rcpp::log(probk_vec);
    NumericVector xi(ploidy + 1);
    for (int k = 0; k <= ploidy; k++) {
      xi(k) = xi_double((double)k / (double)ploidy, seq, bias);
    }

    double sumi; // denominator to get wik for each i.
    NumericVector wvec(ploidy + 2);
    for (int i = 0; i < nind; i++) {
      for (int k = 0; k <= ploidy; k++) {
        wvec(k) = log(1.0 - out_prop) + lprobk_vec(k) + dbetabinom_double(refvec(i), sizevec(i), xi(k), od, true);
      }
      wvec(ploidy + 1) = log(out_prop) + doutdist(refvec(i), sizevec(i), true);
      sumi = log_sum_exp(wvec);
      wvec = Rcpp::exp(wvec - sumi);
      wik_mat(i, _) = wvec;
    }
  }

  return wik_mat;
}




//' Log-likelihood that \code{\link{flexdog}} maximizes when
//' outliers are present.
//'
//' @inheritParams flexdog_full
//' @param probk_vec The kth element is the prior probability of genotype k (when starting to count from 0).
//' @param out_prop The probability of being an outlier.
//'
//' @return A double. The \code{\link{flexdog}} objective when
//'     \code{outliers = TRUE}.
//'
//' @author David Gerard
//'
//' @seealso \code{\link{flexdog_obj}} for the objective function without outliers.
//'
// [[Rcpp::export]]
double flexdog_obj_out(NumericVector probk_vec,
                       double out_prop,
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
    Rcpp::stop("flexdog_obj_out: sizevec and refvec must have the same length.");
  }
  if (probk_vec.length() != ploidy + 1) {
    Rcpp::stop("flexdog_obj_out: probk_vec must have length ploidy + 1.");
  }

  // Calculate the posterior probability of each genotype -------------------
  NumericVector lprobk_vec = Rcpp::log(probk_vec);
  NumericVector xi(ploidy + 1);
  for (int k = 0; k <= ploidy; k++) {
    xi(k) = xi_double((double)k / (double)ploidy, seq, bias);
  }

  // Calculate likelihood ---------------------------------------------------
  double obj = 0.0;
  NumericVector wvec(ploidy + 2);
  if (out_prop < TOL) { // only non-outliers contribute
    for (int i = 0; i < nind; i++) {
      for (int k = 0; k <= ploidy; k++) {
        wvec(k) = lprobk_vec(k) + dbetabinom_double(refvec(i), sizevec(i), xi(k), od, true);
      }
      wvec(ploidy + 1) = R_NegInf;
      obj = obj + log_sum_exp(wvec);
    }
  } else if (1.0 - out_prop < TOL) { // only outliers contribute
    double out_val = 0.0;
    for (int i = 0; i < nind; i++) {
      out_val = doutdist(refvec(i), sizevec(i), true);
      obj = obj + out_val;
    }
  } else {
    for (int i = 0; i < nind; i++) {
      for (int k = 0; k <= ploidy; k++) {
        wvec(k) = log(1.0 - out_prop) + lprobk_vec(k) + dbetabinom_double(refvec(i), sizevec(i), xi(k), od, true);
      }
      wvec(ploidy + 1) = log(out_prop) + doutdist(refvec(i), sizevec(i), true);
      obj = obj + log_sum_exp(wvec);
    }
  }

  // Penalties --------------------------------------------------------------
  obj = obj + pen_bias(bias, mean_bias, var_bias);
  obj = obj + pen_seq_error(seq, mean_seq, var_seq);
  obj = obj + pen_seq_error(od, mean_od, var_od);
  return obj;
}

