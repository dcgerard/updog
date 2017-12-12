#include "mupdog.h"

//' Adjusts allele dosage \code{p} by the sequencing error rate \code{eps}.
//'
//' @param p The allele dosage.
//' @param eps The sequencing error rate.
//'
//' @author David Gerard
// [[Rcpp::export]]
double eta_double(double p, double eps) {
  if (p < -TOL or p > 1.0 + TOL) {
    Rcpp::stop("p must be between 0 and 1");
  }
  if (eps < -TOL or eps > 1.0 + TOL) {
    Rcpp::stop("eps must be between 0 and 1");
  }

  return p * (1 - eps) + (1 - p) * eps;
}

//' Adjusts allele dosage \code{p} by the sequencing error rate \code{eps}.
//'
//' @inheritParams xi_fun
//'
//' @author David Gerard
// [[Rcpp::export]]
NumericVector eta_fun(NumericVector p, double eps) {
  int n = p.length();
  NumericVector eta(n);

  for (int i = 0; i < n; i++) {
    eta(i) = eta_double(p(i), eps);
  }

  return eta;
}

//' Adjusts allele dosage \code{p} by the sequencing error rate \code{eps} and the allele bias \code{h}.
//'
//' @param p A vector of allele dosages.
//' @param eps The sequencing error rate. Must be of length 1.
//' @param h The allele bias. Must be of length 1.
//'
//' @author David Gerard
// [[Rcpp::export]]
NumericVector xi_fun(NumericVector p, double eps, double h) {
  if (h < -TOL) {
    Rcpp::stop("h must be greater than or equal to 0.");
  }

  NumericVector eta_vec = eta_fun(p, eps);

  int n = p.length();
  NumericVector xi(n);

  for (int i = 0; i < n; i++) {
    xi(i) = eta_vec(i) / (h * (1 - eta_vec(i)) + eta_vec(i));
  }

  return xi;
}



