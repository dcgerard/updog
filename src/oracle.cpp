#include "mupdog.h"

// functions for oracle misclassification rates

//' Calculate oracle misclassification error.
//'
//' Given knowledge of the parameters, we calculate the expected misclassification error,
//' where the expectation is taken over both the data generation process and the allele-distribution.
//' This is an ideal level of misclassification error.
//'
//' To come up with \code{dist}, you need some additional assumptions.
//' For example, if the population is in Hardy-Weinberg equilibrium and
//' the allele frequency is \code{alpha} then you could calculate
//' \code{dist} using the R code: \code{dbinom(x = 0:ploidy, size = ploidy, prob = alpha)}.
//' Alternatively, if you know the genotypes of the individual's two parents are, say,
//' \code{ref_count1} and \code{ref_count2}, then you could use the \code{\link[updog]{get_q_array}}
//' function from the updog package: \code{get_q_array(ploidy)[ref_count1 + 1, ref_count2 + 1, ]}.
//'
//' @param n The read-depth.
//' @param ploidy The ploidy of the individual.
//' @param seq The sequencing error rate.
//' @param bias The allele-bias
//' @param od The overdispersion parameter.
//' @param dist The distribution of the alleles.
//'
//'
//' @author David Gerard
//'
//' @export
//'
// [[Rcpp::export]]
double oracle_mis(int n, int ploidy, double seq, double bias, double od, NumericVector dist) {

  if (std::abs(Rcpp::sum(dist) - 1) > TOL) {
    Rcpp::stop("oracle_miss: elements in dist must sum to one.");
  }
  if (dist.length() != ploidy + 1) {
    Rcpp::stop("oracle_miss: dist must have the same length as ploidy + 1.");
  }

  NumericVector ldist = Rcpp::log(dist);

  // Calculate mean values
  NumericVector xi_vec(ploidy + 1);
  for (int i = 0; i <= ploidy; i++) {
    xi_vec(i) = xi_double((double)i / (double)ploidy, seq, bias);
  }

  // Find x's that define the boundaries of the bins
  double val1 = R_NegInf;
  double val2;
  int bin_of_0 = -1; // bin that x = 0 belongs to

  for (int k = 0; k <= ploidy; k++) {
    val2 = dbetabinom_double(0, n, xi_vec(k), od, true) + ldist(k);
    if (val2 > val1) {
      bin_of_0 = k;
      val1 = val2;
    }
  }

  if (bin_of_0 == -1) {
    Rcpp::stop("oracle_miss: bin_of_0 not initialized.");
  }

  // Sum up the max probabilities.
  int current_bin = bin_of_0;
  NumericVector lpbin(ploidy + 1, R_NegInf);

  for (int x = 0; x <= n; x++) {
    if (current_bin < ploidy) {
      val1 = dbetabinom_double(x, n, xi_vec(current_bin), od, true) + ldist(current_bin);
      val2 = dbetabinom_double(x, n, xi_vec(current_bin + 1), od, true) + ldist(current_bin + 1);
      if (val2 > val1) {
        current_bin++;
        val1 = val2;
      }
    } else {
      val1 = dbetabinom_double(x, n, xi_vec(current_bin), od, true) + ldist(current_bin);
    }
    lpbin(current_bin) = log_sum_exp_2(lpbin(current_bin), val1);
  }

  double miss_error = 1.0 - std::exp(log_sum_exp(lpbin));

  return miss_error;
}
