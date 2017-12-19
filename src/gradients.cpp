#include "mupdog.h"

//' Gradient for \code{\link{obj_for_mu_sigma2}} with respect for \code{mu} and \code{sigma2}.
//'
//' @inheritParams obj_for_mu_sigma2
//'
//' @author David Gerard
// [[Rcpp::export]]
NumericVector grad_for_mu_sigma2(arma::Col<double> mu, arma::Col<double> sigma2, NumericMatrix phifk_mat,
                                 arma::Mat<double> cor_inv, NumericMatrix log_bb_dense) {
  // Check input --------------------------------------------------------------------------------
  int nind = log_bb_dense.nrow();
  int ploidy = log_bb_dense.ncol() - 1;
  if (mu.n_elem != (unsigned int)nind) {
    Rcpp::stop("obj_for_mu_sigma2: mu needs to have the same length as the number of columns in log_bb_dense.");
  }
  if (sigma2.n_elem != (unsigned int)nind) {
    Rcpp::stop("obj_for_mu_sigma2: sigma2 needs to have the same length as the number of columns in log_bb_dense.");
  }
  if (phifk_mat.nrow() != nind) {
    Rcpp::stop("obj_for_mu_sigma2: phifk_mat and log_bb_dense needs to have the same dimensions.");
  }
  if (phifk_mat.ncol() != ploidy + 2) {
    Rcpp::stop("obj_for_mu_sigma2: phifk_mat needs to have ploidy+2 columns.");
  }
  if (cor_inv.n_cols != (unsigned int)nind) {
    Rcpp::stop("obj_for_mu_sigma2: cor_inv needs to have the same number of rows as log_bb_dense.");
  }
  if (cor_inv.n_rows != (unsigned int)nind) {
    Rcpp::stop("obj_for_mu_sigma2: cor_inv needs to have the same number of columns as rows.");
  }

  NumericVector grad(2 * nind); // first n for mu, next nf or sigma2

  // contribution by likelihood ------------------------------------------------
  double current_weight = 0.0;
  double current_sigma_weight = 0.0;
  double q1 = 0.0;
  double q2 = 0.0;
  double d1 = 0.0;
  double d2 = 0.0;
  for (int i = 0; i < nind; i++) {
    for (int k = 0; k <= ploidy; k++) {
      if (!R_IsNA(log_bb_dense(i, k))) {
        q1 = (phifk_mat(i, k) - mu(i)) / std::sqrt(sigma2(i));
        q2 = (phifk_mat(i, k + 1) - mu(i)) / std::sqrt(sigma2(i));
        d1 = R::dnorm4(q1, 0.0, 1.0, false);
        d2 = R::dnorm4(q2, 0.0, 1.0, false);
        current_weight = d1 / std::sqrt(sigma2(i)) - d2 / std::sqrt(sigma2(i));
        grad(i) = grad(i) + current_weight * log_bb_dense(i, k);

        if ((q1 != R_PosInf) & (q1 != R_NegInf) & (q2 != R_PosInf) & (q2 != R_NegInf)) {
          current_sigma_weight = d1 * q1 / (2.0 * sigma2(i)) - d2 * q2 / (2.0 * sigma2(i));
        } else if ((q1 == R_NegInf) & (q2 != R_NegInf) & (q2 != R_PosInf)) {
          current_sigma_weight = -d2 * q2 / (2.0 * sigma2(i));
        } else if ((q1 != R_NegInf) & (q1 != R_PosInf) & (q2 == R_PosInf)){
          current_sigma_weight = d1 * q1 / (2.0 * sigma2(i));
        } else {
          current_sigma_weight = 0.0;
        }
        grad(i + nind) = grad(i + nind) + current_sigma_weight * log_bb_dense(i, k);
      }
    }
  }

  // contribution by prior ----------------------------------------------------
  arma::Col<double> sig2_pen = arma::diagvec(cor_inv) / -2.0 + 1.0 / (2.0 * sigma2);
  arma::Col<double> mu_pen = -1.0 * cor_inv * mu;
  for (int i = 0; i < nind; i++) {
    grad(i) = grad(i) + mu_pen(i);
    grad(i + nind) = grad(i + nind) + sig2_pen(i);
  }

  return grad;
}


//' Gradient for \code{\link{obj_for_mu_sigma2_wrapper}} with respect for \code{muSigma2}
//' and a wrapper for \code{\link{grad_for_mu_sigma2}}
//'
//' @inheritParams obj_for_mu_sigma2
//' @param muSigma2 A vector. The first half are mu and the second half are sigma2.
//'
//' @author David Gerard
// [[Rcpp::export]]
NumericVector grad_for_mu_sigma2_wrapper(arma::Col<double> muSigma2, NumericMatrix phifk_mat,
                                         arma::Mat<double> cor_inv, NumericMatrix log_bb_dense) {
  int nind = muSigma2.n_elem / 2;
  NumericVector grad = grad_for_mu_sigma2(muSigma2.head_rows(nind), muSigma2.tail_rows(nind),
                                          phifk_mat, cor_inv, log_bb_dense);
  return grad;
}


