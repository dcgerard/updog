#include <Rcpp.h>
using namespace Rcpp;

//' Adjusts allele dosage \code{p} by the sequencing error rate \code{eps}.
//'
//' @param p The allele dosage.
//' @param eps The sequencing error rate.
// [[Rcpp::export]]
NumericVector eta_fun(double p, double eps) {
  return p * (1 - eps) + (1 - p) * eps;
}
