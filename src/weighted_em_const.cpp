#include "mupdog.h"

// Functions for solving weighted EM with a uniform mixing component.


//' Objective function optimized by \code{\link{uni_em_const}}.
//'
//' @inheritParams uni_em_const
//' @param pivec The current parameters.
//'
//' @author David Gerard
//'
//' @return The objective optimized by \code{\link{uni_em_const}} during
//'     that separate unimodal EM algorithm.
//'
// [[Rcpp::export]]
double uni_obj_const(arma::vec pivec,
                     double alpha,
                     arma::vec weight_vec,
                     arma::mat lmat,
                     long double lambda) {
  arma::vec lpi = (1.0 - alpha) * lmat.t() * pivec + (alpha / (double)weight_vec.n_elem);
  double obj = 0.0;
  for (int k = 0; k < weight_vec.n_elem; k++) {
    if ((weight_vec(k) > TOL) && (lpi(k) > TOL)) {
      obj = obj + weight_vec(k) * std::log(lpi(k));
    } else if ((weight_vec(k) > TOL * 1000) && (lpi(k) < TOL)) {
      obj = R_NegInf;
      break;
    } else {
      // do nothing.
    }
  }

  // Add the penalty ----
  double pen = 0.0;
  if (lambda > TOL) {
    pen = lambda * arma::sum(arma::log(pivec));
  }

  return obj + pen;
}

//' EM algorithm to fit weighted ash objective with a uniform
//' mixing component.
//'
//' Solves the following optimization problem
//' \deqn{\max_{\pi} \sum_k w_k \log(\alpha / (K+1) + (1 - \alpha)\sum_j \pi_j \ell_jk).}
//' It does this using a weighted EM algorithm.
//'
//' @param weight_vec A vector of weights. Each element of \code{weight_vec} corresponds
//'     to a column of \code{lmat}.
//' @param lmat A matrix of inner weights. The columns are the "individuals" and the rows are the "classes."
//' @param pi_init The initial values of \code{pivec}. Each element of \code{pi_init}
//'     corresponds to a row of \code{lmat}.
//' @param alpha The mixing weight for the uniform component.
//'     This should be small (say, less tahn 10^-3).
//' @param itermax The maximum number of EM iterations to take.
//' @param obj_tol The objective stopping criterion.
//' @param lambda The penalty on the pi's. Should be greater than 0 and really really small.
//'
//'
//' @return A vector of numerics. The update of \code{pivec} in
//'     \code{\link{flexdog_full}}.
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
arma::vec uni_em_const(arma::vec weight_vec,
                       arma::mat lmat,
                       arma::vec pi_init,
                       double alpha,
                       long double lambda,
                       int itermax,
                       double obj_tol) {
  // check input ----------------------------------------------
  int nind = weight_vec.n_elem;
  int nclass = pi_init.n_elem;
  if (obj_tol < TOL) {
    Rcpp::stop("uni_em_const: obj_tol should be greater than 0.");
  }
  if (itermax < 0) {
    Rcpp::stop("uni_em_const: itermax should be greater than or equal to 0.");
  }
  if (lmat.n_rows != nclass) {
    Rcpp::stop("uni_em_const: lmat should have pi_init.n_elem rows.");
  }
  if (lmat.n_cols != nind) {
    Rcpp::stop("uni_em_const: lmat should have weight_vec.n_elem columns.");
  }
  if (lambda < 0.0) {
    Rcpp::stop("uni_em_const: lambda cannot be negative.");
  }
  if ((alpha < 0.0) | (alpha > 1.0 - TOL)) {
    Rcpp::stop("uni_em_const: alpha should be in [0, 1).");
  }

  // Run EM ---------------------------------------------------
  int index       = 0;
  double err      = obj_tol + 1.0;
  arma::vec pivec = pi_init;
  double obj      = uni_obj_const(pivec, alpha, weight_vec, lmat, lambda);
  double old_obj  = obj;
  double lsum     = 0.0;
  arma::mat etamat(nclass, nind);
  arma::vec nvec(nclass);

  while ((index < itermax) & (err > obj_tol)) {
    old_obj = obj;
    // get eta_jk -----------------
    for (int k = 0; k < nind; k++) {
      lsum = 0.0;
      for (int j = 0; j < nclass; j++) {
        etamat(j, k) = (1.0 - alpha) * pivec(j) * lmat(j, k);
        lsum = lsum + etamat(j, k);
      }
      lsum = lsum + alpha / (double)nind;
      for (int j = 0; j < nclass; j++) {
        etamat(j, k) = etamat(j, k) / lsum;
      }
    }
    // get n_j's --------------------------------
    nvec = etamat * weight_vec + lambda;
    // normalize to get pi_j's ------------------
    pivec = nvec / arma::sum(nvec);
    // calculate objective and update stopping criteria
    obj = uni_obj_const(pivec, alpha, weight_vec, lmat, lambda);
    if (obj < old_obj - TOL) {
      Rcpp::stop("uni_em: Objective is not increasing.\n");
    }
    err = std::abs(obj - old_obj);
    index++;
  }
  return pivec;
}

//' Convolution between two discrete probability mass functions
//' with support on 0:K.
//'
//' @author David Gerard
//'
//' @param x The first probability vector. The ith element is the
//'     probability of i - 1.
//' @param y The second probability vector. The ith element is the
//'     probability of i - 1.
//'
//' @return A vector that is the convolution of \code{x} and
//'     \code{y}. The ith element is the probability of i - 1.
//'
//' @export
// [[Rcpp::export]]
arma::vec convolve(arma::vec x, arma::vec y) {
  if (x.n_elem != y.n_elem) {
    Rcpp::stop("convolve: x and y should have the same number of values.");
  }

  int nval = x.n_elem;
  arma::mat joint(nval, nval);

  for (int i = 0; i < nval; i++) {
    for (int j = 0; j < nval; j++) {
      joint(i, nval - j - 1) = x(i) * y(j);
    }
  }

  arma::vec conv(2 * nval - 1);
  for (int i = 0; i < (2 * nval - 1); i++) {
    conv(2 * nval - 2 - i) = arma::accu(joint.diag(i - nval + 1));
  }

  return(conv);
}



//' Objective function when doing Brent's method in
//' \code{\link{update_pp}} when one parent only has
//' two mixing components.
//'
//' @param firstmixweight The mixing weight of the first component.
//' @param probmat The rows index the components and the columns
//'     index the segregation amount. Should only have two rows.
//' @param pvec The distribution of the other parent.
//' @param weight_vec The weights for each element.
//' @param alpha The mixing weight on the uniform component.
//'
//' @return The objective value, as calculated by taking a
//'     convolution using \code{\link{convolve}} of the mixing
//'     distribution and \code{pvec}, then putting that
//'     probability distribution through \code{\link{f1_obj}}.
//'
//' @author David Gerard
// [[Rcpp::export]]
double pp_brent_obj(double firstmixweight,
                    arma::mat probmat,
                    arma::vec pvec,
                    arma::vec weight_vec,
                    double alpha) {
  int ploidy = weight_vec.n_elem - 1;

  if (probmat.n_rows != 2) {
    Rcpp::stop("pp_brent_obj: probmat should have two rows.");
  }
  if (probmat.n_cols != (ploidy / 2 + 1)) {
    Rcpp::stop("pp_brent_obj: probmat should have ploidy / 2 + 1 columns.");
  }
  if (probmat.n_cols != pvec.n_elem) {
    Rcpp::stop("pp_brent_obj: probmat.n_cols should equal pvec.n_elem.");
  }
  if ((alpha < 0.0) | (alpha > 1.0 - TOL)) {
    Rcpp::stop("pp_brent_obj: alpha should be in [0, 1)");
  }
  if ((firstmixweight < 0.0) | (firstmixweight > 1.0)) {
    Rcpp::stop("pp_brent_obj: firstmixweight should be in [0, 1]");
  }

  arma::vec pvec_new = firstmixweight * probmat.row(0).t() +
    (1.0 - firstmixweight) * probmat.row(1).t();

  arma::vec pvec_final = convolve(pvec_new, pvec);

  double obj = f1_obj(alpha, pvec_final, weight_vec);

  return obj;
}
