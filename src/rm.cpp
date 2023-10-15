// code for real_to_simplex() and dreal_to_simplex()
#include "mupdog.h"

// from util_funs.cpp
double expit(double x);

//' Convert real line to simplex using Stan technique
//'
//' @param y A vector of numbers of length K-1
//'
//' @return A vector on the simplex of length K
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
arma::vec real_to_simplex(const arma::vec y) {
  int K = y.n_elem + 1;
  arma::vec x(K);

  double recsum = 0.0;
  for (int k = 0; k < K - 1; k++) {
    x[k] = (1.0 - recsum) * expit(y[k] + std::log(1.0 / ((double)K - ((double)k + 1))));
    recsum += x[k];
  }

  x[K - 1] = 1.0 - recsum;

  return x;
}

//' Objective function for random mating m step in EM
//'
//' @param y The real-valued parameterization of the simplex
//' @param weight_vec The current weights
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
double obj_rm(arma::vec y, arma::vec weight_vec) {
  arma::vec p = real_to_simplex(y);
  arma::vec q = arma::conv(p, p);
  arma::vec lvec = weight_vec % arma::log(q);
  double obj = arma::sum(lvec);
  return obj;
}

//' Derivative of \code{\link{real_to_simplex}()} with respect to \code{y}.
//'
//' @param y A numeric matrix.
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
arma::mat dreal_to_simplex_dy(const arma::vec y) {
  int K = y.n_elem + 1;

  arma::mat jacob(K, K - 1);
  double recsum;
  double drecsum;
  double zk;
  for (int j = 0; j < K - 1; j++) {
    recsum = 0.0;
    drecsum = 0.0;
    for (int i = 0; i < K; i++) {
      if (i < K - 1) {
        zk = expit(y[i] + std::log(1.0 / ((double)K - ((double)i + 1))));
      } else{
        zk = expit(std::log(1.0 / ((double)K - ((double)i + 1))));
      }
      if (i < j) {
        jacob(i, j) = 0.0;
      } else if (i == j) {
        jacob(i, j) = (1 - recsum) * zk * (1.0 - zk);
      } else {
        jacob(i, j) = - drecsum * zk;
      }
      recsum += (1.0 - recsum) * zk;
      drecsum += jacob(i, j);
    }
  }

  return jacob;
}


//' Derivative of q = \code{stats::convolve(p, rev(p), type = "open")}
//' with respect to p.
//'
//' @param p The gamete frequencies
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
arma::mat dq_dp(const arma::vec p) {
  int khalf = p.n_elem; // ploidy / 2 + 1
  int K = 2 * (khalf - 1); // ploidy
  arma::mat jacob(K + 1, khalf);
  for (int k = 0; k <= K; k++) {
    for (int j = 0; j < khalf; j++) {
      if (j >= std::max(0, k - K/2) && j <= std::min(k, K/2)) {
        jacob(k, j) = 2.0 * p(k - j);
      }
    }
  }
  return jacob;
}

//' Derivative of \code{\link{rm_llike}()} with respect to q.
//'
//' @param q The genotype frequencies
//'
//' @noRd
// [[Rcpp::export]]
arma::vec drmlike_dq(arma::vec q, arma::vec weight_vec) {
  arma::vec jacob = weight_vec / q;
  return jacob;
}

//' Derivative of obj_rm() with respect to y
//'
//' @inheritParams obj_rm
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
arma::vec dobjrm_dy(arma::vec y, arma::vec weight_vec) {
  // intermediate values ----
  arma::vec p = real_to_simplex(y);
  arma::vec q = arma::conv(p, p);

  // intermediate derivatives ----
  arma::mat j_dp_dy = dreal_to_simplex_dy(y);
  arma::mat j_dq_dp = dq_dp(p);
  arma::vec j_df_dq = drmlike_dq(q, weight_vec);

  // gradient ----
  arma::vec jacob = j_dp_dy.t() * j_dq_dp.t() * j_df_dq;

  Rcpp::Rcout << jacob << std::endl;

  return jacob;
}



