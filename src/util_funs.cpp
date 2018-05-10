#include "mupdog.h"

//' Adjusts allele dosage \code{p} by the sequencing error rate \code{eps}.
//'
//' @param p The allele dosage.
//' @param eps The sequencing error rate.
//'
//' @return The probability of a reference reed adjusted by the
//'     sequencing error rate.
//'
//' @author David Gerard
// [[Rcpp::export]]
double eta_double(double p, double eps) {
  if ((p < -TOL) or (1.0 - p < -TOL)) {
    Rcpp::stop("eta_double: p must be between 0 and 1");
  }
  if ((eps < -TOL) or (1.0 - eps < -TOL)) {
    Rcpp::stop("eta_double: eps must be between 0 and 1");
  }

  return p * (1 - eps) + (1 - p) * eps;
}

//' Adjusts allele dosage \code{p} by the sequencing error rate \code{eps}.
//'
//' @inheritParams xi_fun
//'
//' @return A vector of probabilities of a reference read adjusted
//'     by the sequencing error rate.
//'
//' @author David Gerard
// [[Rcpp::export]]
NumericVector eta_fun(NumericVector p, NumericVector eps) {
  int n = p.length();
  if ((eps.length() != n) & (eps.length() != 1)) {
    Rcpp::stop("eta_fun: eps must either have length 1 or be the same length as p.");
  }

  NumericVector eta(n);
  double eps_current;

  for (int i = 0; i < n; i++) {
    if (eps.length() == n) {
      eps_current = eps(i);
    }
    else {
      eps_current = eps(0);
    }

    eta(i) = eta_double(p(i), eps_current);
  }

  return eta;
}



//' Adjusts allele dosage \code{p} by the sequencing error rate \code{eps} and the allele bias \code{h}.
//'
//' @param p The allele dosage.
//' @param eps The sequencing error rate.
//' @param h The allele bias.
//'
//' @return The probability of a reference read adjusted by both the allele
//'     bias and the sequencing error rate.
//'
//' @author David Gerard
// [[Rcpp::export]]
double xi_double(double p, double eps, double h) {
  if (h < -TOL) {
    Rcpp::stop("xi_double: h must be greater than or equal to 0.");
  }
  double eta = eta_double(p, eps);
  double xi  = eta / (h * (1.0 - eta) + eta);
  return xi;
}

//' Adjusts allele dosage \code{p} by the sequencing error rate \code{eps} and the allele bias \code{h}.
//'
//' @param p A vector of allele dosages.
//' @param eps The sequencing error rate. Must either be of length 1
//'     or the same length as p.
//' @param h The allele bias. Must either be of length 1 or the same length
//'     as p.
//'
//' @return A vector of prababilities of the reference read adjusted
//'     by both the sequencing error rate and the allele bias.
//'
//' @author David Gerard
// [[Rcpp::export]]
NumericVector xi_fun(NumericVector p, NumericVector eps, NumericVector h) {
  int n = p.length();
  if ((eps.length() != n) & (eps.length() != 1)) {
    Rcpp::stop("xi_fun: eps must either have length 1 or the same length as x.");
  }
  if ((h.length() != n) & (h.length() != 1)) {
    Rcpp::stop("xi_fun: h must either have length 1 or the same length as x.");
  }

  NumericVector xi(n);
  double eps_current;
  double h_current;

  for (int i = 0; i < n; i++) {
    if (eps.length() == n) {
      eps_current = eps(i);
    }
    else {
      eps_current = eps(0);
    }

    if (h.length() == n) {
      h_current = h(i);
    }
    else {
      h_current = h(0);
    }

    xi(i) = xi_double(p(i), eps_current, h_current);
  }

  return xi;
}


//' Log-sum-exponential trick.
//'
//' @param x A vector to log-sum-exp.
//'
//' @return The log of the sum of the exponential
//'     of the elements in \code{x}.
//'
//' @author David Gerard
// [[Rcpp::export]]
double log_sum_exp(NumericVector x) {
  double max_x = Rcpp::max(x);
  double lse; // the log-sum-exp
  if (max_x == R_NegInf) { // if all -Inf, need to treat this special to avoid -Inf + Inf.
    lse = R_NegInf;
  }
  else {
    lse = max_x + std::log(Rcpp::sum(Rcpp::exp(x - max_x)));
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
// [[Rcpp::export]]
double log_sum_exp_2(double x, double y) {
  double z = std::max(x, y);
  double finalval;
  if (z == R_NegInf) {
    finalval = R_NegInf;
  } else {
    finalval = std::log(std::exp(x - z) + std::exp(y - z)) + z;
  }
  return finalval;
}

//' The logit function.
//'
//' @param x A double between 0 and 1.
//'
//' @return The logit of \code{x}.
//'
//' @author David Gerard
// [[Rcpp::export]]
double logit(double x) {
  if ((x < TOL) | ((1.0 - x) < TOL)) {
    Rcpp::stop("logit: x must be between 0 and 1.");
  }
  double lv = std::log(x / (1.0 - x));
  return lv;
}

//' The expit (logistic) function.
//'
//' @param x A double.
//'
//' @return The expit (logistic) of \code{x}.
//'
//' @author David Gerard
// [[Rcpp::export]]
double expit(double x) {
  double ev = 1.0 / (1.0 + std::exp(-x));
  return ev;
}
