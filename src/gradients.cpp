#include "mupdog.h"

//' Gradient for \code{\link{obj_for_mu_sigma2}} with respect for \code{mu} and \code{sigma2}.
//'
//' @inheritParams obj_for_mu_sigma2
//'
//' @return A vector of length 2 * nind of numerics.
//'     The first element n elements are the partial derivatives with respect
//'     to \code{mu} and the second n elements are the
//'     partial derivatives with respect to \code{sigma2}
//'     in \code{\link{obj_for_mu_sigma2}}.
//'
//' @author David Gerard
//'
//' @keywords internal
//' @noRd
//'
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

        if ((q1 != R_PosInf) && (q1 != R_NegInf) && (q2 != R_PosInf) && (q2 != R_NegInf)) {
          current_sigma_weight = d1 * q1 / (2.0 * sigma2(i)) - d2 * q2 / (2.0 * sigma2(i));
        } else if ((q1 == R_NegInf) && (q2 != R_NegInf) && (q2 != R_PosInf)) {
          current_sigma_weight = -d2 * q2 / (2.0 * sigma2(i));
        } else if ((q1 != R_NegInf) && (q1 != R_PosInf) && (q2 == R_PosInf)){
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
//'
//' @inherit grad_for_mu_sigma2 return
//' @inheritParams obj_for_mu_sigma2
//' @param muSigma2 A vector. The first half are mu and the second half are sigma2.
//'
//' @keywords internal
//' @noRd
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


//' Derivative of \deqn{-log(h) - (log(h) - \mu_h)^2 / (2\sigma_h^2)} with respect
//' to \eqn{h}.
//'
//' @param h The current bias parameter.
//' @param mu_h The mean of the log-bias.
//' @param sigma2_h The variance of the log-bias.
//'
//' @seealso \code{\link{pen_bias}} which this is a derivative for.
//'
//' @return A double.
//'
//' @keywords internal
//' @noRd
//'
//' @author David Gerard
// [[Rcpp::export]]
double dpen_dh(double h, double mu_h, double sigma2_h) {
  double deriv;
  if (arma::is_finite(sigma2_h)) {
    deriv = -1.0 * (1.0 + (std::log(h) - mu_h) / sigma2_h) / h;
  } else {
    deriv = 0.0;
  }
  return deriv;
}


//' Derivative of \deqn{-log(\epsilon(1 - \epsilon)) - (logit(\epsilon) - \mu_{\epsilon})^2 / (2\sigma_{\epsilon}^2)}
//' with respect to \eqn{\epsilon}.
//'
//' @param eps The current sequencing error rate.
//' @param mu_eps The mean of the logit of the sequencing error rate.
//' @param sigma2_eps The variance of the logit of the sequencing error rate.
//'
//' @seealso \code{\link{pen_seq_error}} which this is a derivative for.
//'
//' @return A double.
//'
//' @keywords internal
//' @noRd
//'
//' @author David Gerard
// [[Rcpp::export]]
double dpen_deps(double eps, double mu_eps, double sigma2_eps) {
  double deriv;

  if (arma::is_finite(sigma2_eps)) {
    deriv = -1.0 * (1.0 - 2.0 * eps + (logit(eps) - mu_eps) / sigma2_eps) / (eps * (1.0 - eps));
  } else {
    deriv = 0.0;
  }
  return deriv;
}


//' Derivative of the log-beta density with
//' respect to c where \eqn{c = (1 - \tau)/\tau}
//' where \eqn{\tau} is the overdispersion parameter.
//'
//' @param x The number of successes observed
//' @param n The total number of trials observed.
//' @param xi The mean of the beta-binomial.
//' @param c \eqn{(1 - \tau)/\tau} where \eqn{\tau}
//'     is the overdispersion parameter.
//'
//' @return A double.
//'
//' @keywords internal
//' @noRd
//'
//' @seealso \code{\link{dbetabinom_double}}, \code{\link{dlbeta_dtau}},
//'     \code{\link{dc_dtau}}.
//' @author David Gerard
// [[Rcpp::export]]
double dlbeta_dc(int x, int n, double xi, double c) {
  double deriv = -1.0 * xi * R::digamma(xi * c) -
    (1.0 - xi) * R::digamma((1.0 - xi) * c) +
    R::digamma(c) +
    xi * R::digamma((double)x + xi * c) +
    (1.0 - xi) * R::digamma((double)n - (double)x + (1.0 - xi) * c) -
    R::digamma((double)n + c);

  return(deriv);
}


//' Derivative of \eqn{c = (1 - \tau) / \tau} with respect to \eqn{\tau}.
//'
//' @param tau The overdispersion parameter.
//'
//' @return A double.
//'
//' @keywords internal
//' @noRd
//'
//' @seealso \code{\link{dlbeta_dc}}, \code{\link{dlbeta_dtau}}
//'
//' @author David Gerard
// [[Rcpp::export]]
double dc_dtau(double tau) {
  double deriv = -1.0 / std::pow(tau, 2.0);
  return deriv;
}

//' Derivative of the log-beta-binomial density with respect to the
//' overdispersion parameter.
//'
//' @param x The number of successes.
//' @param n The number of trials.
//' @param tau The overdispersion parameter.
//' @param p The allele dosage.
//' @param eps The sequencing error rate
//' @param h The bias parameter.
//'
//' @return A double.
//'
//' @keywords internal
//' @noRd
//'
//' @seealso \code{\link{dlbeta_dc}}, \code{\link{dc_dtau}},
//'     \code{\link{dbetabinom_double}}.
//'
//' @author David Gerard
// [[Rcpp::export]]
double dlbeta_dtau(int x, int n, double p, double eps, double h, double tau) {
  double xi = xi_double(p, eps, h);
  double dlbetadc = dlbeta_dc(x, n, xi, (1.0 - tau) / tau);
  double dcdtau   = dc_dtau(tau);
  double deriv = dlbetadc * dcdtau;
  return deriv;
}


//' Derivative of the log-betabinomial density with respect to the
//' mean of the underlying beta.
//'
//' @return A double.
//'
//' @inheritParams dlbeta_dtau
//' @param xi The mean of the underlying beta.
//'
//' @keywords internal
//' @noRd
//'
//' @author David Gerard
// [[Rcpp::export]]
double dlbeta_dxi(int x, int n, double xi, double tau) {
  double c = (1.0 - tau) / tau;
  double deriv = c * R::digamma((double)x + xi * c) -
    c * R::digamma((double)n - (double)x + (1.0 - xi) * c) -
    c * R::digamma(xi * c) +
    c * R::digamma((1.0 - xi) * c);
  return deriv;
}

//' Derivative of xi-function with respect to bias parameter.
//'
//' @return A double.
//'
//' @param p The dosage (between 0 and 1).
//' @param eps The sequencing error rate.
//' @param h The bias parameter.
//'
//' @keywords internal
//' @noRd
//'
//' @author David Gerard
// [[Rcpp::export]]
double dxi_dh(double p, double eps, double h) {
  double f = p * (1.0 - eps) + eps * (1.0 - p);
  double deriv = -1.0 * f * (1.0 - f) / std::pow(h * (1.0 - f) + f, 2.0);
  return deriv;
}

//' Derivative of log-betabinomial density with respect to bias parameter.
//'
//' @inheritParams dlbeta_dtau
//'
//' @return A double.
//'
//' @keywords internal
//' @noRd
//'
//' @author David Gerard
// [[Rcpp::export]]
double dlbeta_dh(int x, int n, double p, double eps, double h, double tau) {
  double xi = xi_double(p, eps, h);
  double dlbetadxi = dlbeta_dxi(x, n, xi, tau);
  double dxidh = dxi_dh(p, eps, h);
  double deriv = dlbetadxi * dxidh;
  return deriv;
}

//' Derivative of xi with respect to f.
//'
//' @param h The bias parameter.
//' @param f The post-sequencing error rate adjusted probability of an A.
//'
//' @return A double.
//'
//' @keywords internal
//' @noRd
//'
//' @author David Gerard
// [[Rcpp::export]]
double dxi_df(double h, double f) {
  double deriv = h / std::pow(h * (1.0 - f) + f, 2.0);
  return deriv;
}

//' Derivative of f with respect to eps.
//'
//' @param p The allele dosage.
//' @param eps The sequencing error rate.
//'
//' @return A double.
//'
//' @keywords internal
//' @noRd
//'
//' @author David Gerard
// [[Rcpp::export]]
double df_deps(double p, double eps) {
  double deriv = 1.0 - 2.0 * p;
  return deriv;
}

//' Derivative of the log-beta-binomial density with respect to the
//' sequencing error rate.
//'
//' @inheritParams dlbeta_dtau
//'
//' @return A double.
//'
//' @keywords internal
//' @noRd
//'
//' @author David Gerard
// [[Rcpp::export]]
double dlbeta_deps(int x, int n, double p, double eps, double h, double tau) {
  double f = eps * (1.0 - p) + (1.0 - eps) * p;
  double xi = xi_double(p, eps, h);
  double dlbetadxi = dlbeta_dxi(x, n, xi, tau);
  double dxidf = dxi_df(h, f);
  double dfdeps = df_deps(p, eps);
  double deriv = dlbetadxi * dxidf * dfdeps;
  return deriv;
}


//' Gradient for \code{\link{obj_for_eps}}.
//'
//' @return A double.
//'
//' @keywords internal
//' @noRd
//'
//' @inheritParams obj_for_eps
//'
//' @author David Gerard
// [[Rcpp::export]]
NumericVector grad_for_eps(NumericVector parvec,
                   NumericVector refvec,
                   NumericVector sizevec,
                   int ploidy,
                   double mean_bias,
                   double var_bias,
                   double mean_seq,
                   double var_seq,
                   double mean_od,
                   double var_od,
                   NumericMatrix wmat,
                   bool update_bias = true,
                   bool update_seq = true,
                   bool update_od = true) {
  // check input --------------------------------------------------------
  int nind = refvec.length();
  if (sizevec.length() != nind) {
    Rcpp::Rcout << sizevec.length() << std::endl;
    Rcpp::stop("grad_for_eps: sizevec and refvec must have same length.");
  }
  if (wmat.nrow() != nind) {
    Rcpp::Rcout << wmat.nrow() << std::endl;
    Rcpp::stop("grad_for_eps: wmat must have the same number of rows as the length of refved.");
  }
  if (wmat.ncol() != (ploidy + 1)) {
    Rcpp::Rcout << wmat.ncol() << std::endl;
    Rcpp::stop("grad_for_eps: wmat must have ploidy+1 columns.");
  }


  NumericVector grad(3);

  double eps = parvec(0);
  double h   = parvec(1);
  double tau = parvec(2);

  double p;
  for (int i = 0; i < nind; i++) {
    if (!R_IsNA(refvec(i)) && !R_IsNA(sizevec(i))) {
      for (int k = 0; k <= ploidy; k++) {
        p = (double)k / (double)ploidy;
        grad(0) = grad(0) + wmat(i, k) * dlbeta_deps(refvec(i), sizevec(i), p, eps, h, tau);
        grad(1) = grad(1) + wmat(i, k) * dlbeta_dh(refvec(i), sizevec(i), p, eps, h, tau);
        grad(2) = grad(2) + wmat(i, k) * dlbeta_dtau(refvec(i), sizevec(i), p, eps, h, tau);
      }
    }
  }

  // contribution by penalties
  grad(0) = grad(0) + dpen_deps(eps, mean_seq, var_seq);
  grad(1) = grad(1) + dpen_dh(h, mean_bias, var_bias);
  grad(2) = grad(2) + dpen_deps(tau, mean_od, var_od); // can use dpen_deps() for the od as well!

  // set gradient to zero where don't update parameter values.
  if (!update_seq) {
    grad(0) = 0.0;
  }
  if (!update_bias) {
    grad(1) = 0.0;
  }
  if (!update_od) {
    grad(2) = 0.0;
  }

  return grad;
}


//' Gradient for \code{\link{obj_for_weighted_lbb}}.
//'
//' @inheritParams obj_for_weighted_lbb
//'
//' @return A vector of length 2. The first component is the gradient of the mean of the underlying
//'     beta. The second component is the gradient of the overdispersion parameter of the
//'     underlying beta.
//'
//' @keywords internal
//' @noRd
//'
//' @author David Gerard
// [[Rcpp::export]]
NumericVector grad_for_weighted_lbb(NumericVector parvec,
                                    int ploidy,
                                    NumericVector weight_vec) {
  if (parvec.length() != 2) {
    Rcpp::stop("obj_for_weighted_lbb: parvec not of length 2.");
  }
  if (weight_vec.length() != (ploidy + 1)) {
    Rcpp::stop("obj_for_weighted_lbb: weight_vec not of length ploidy + 1.");
  }

  double mu  = parvec(0);
  double tau = parvec(1);
  NumericVector grad(2);
  for (int i = 0; i <= ploidy; i++) {
    grad(0) = grad(0) + weight_vec(i) * dlbeta_dxi(i, ploidy, mu, tau);
    grad(1) = grad(1) + weight_vec(i) * dlbeta_dc(i, ploidy, mu, (1.0 - tau) / tau) * dc_dtau(tau);
  }

  return grad;
}


//' Gradient for \code{\link{obj_for_weighted_lnorm}}.
//'
//' @inheritParams obj_for_weighted_lnorm
//'
//' @return A vector of length 2. The first term is the derivative with respect to the mean,
//'     the second term is the derivative with respect to the standard deviation (not variance).
//'
//' @author David Gerard
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericVector grad_for_weighted_lnorm(NumericVector parvec,
                                      int ploidy,
                                      NumericVector weight_vec) {
  if (parvec.length() != 2) {
    Rcpp::stop("obj_for_weighted_lbb: parvec not of length 2.");
  }
  if (weight_vec.length() != (ploidy + 1)) {
    Rcpp::stop("obj_for_weighted_lbb: weight_vec not of length ploidy + 1.");
  }

  double mu    = parvec(0);
  double sigma = parvec(1);
  NumericVector grad(2);
  NumericVector pvec(ploidy + 1);
  for (int i = 0; i <= ploidy; i++) {
    pvec(i) = R::dnorm(((double)i - mu) / sigma, 0, 1, true);
  }
  double lsum = log_sum_exp(pvec);
  pvec = Rcpp::exp(pvec - lsum);
  double wsum = Rcpp::sum(weight_vec);

  for (int i = 0; i <= ploidy; i++) {
    grad(0) = grad(0) + weight_vec(i) * ((double)i - mu) - wsum * ((double)i - mu) * pvec(i);
    grad(1) = grad(1) + weight_vec(i) * std::pow(((double)i - mu), 2) - wsum * std::pow(((double)i - mu), 2) * pvec(i);
  }

  grad(0) = grad(0) / std::pow(sigma, 2);
  grad(1) = grad(1) / std::pow(sigma, 3);


  return(grad);
}

