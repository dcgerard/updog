#include "mupdog.h"

// functions for oracle misclassification rates

//' Calculate oracle misclassification error rate.
//'
//' Given perfect knowledge of the data generating parameters,
//' \code{oracle_mis} calculates the misclassification error
//' rate, where the error rate is taken over both the data generation
//' process and the allele-distribution.
//' This is an ideal level of the misclassification error rate and
//' any real method will have a larger rate than this. This is a useful
//' approximation when you have a lot of individuals.
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
//' @param bias The allele-bias.
//' @param od The overdispersion parameter.
//' @param dist The distribution of the alleles.
//'
//' @return A double. The oracle misclassification error rate.
//'
//' @author David Gerard
//'
//' @references Gerard, David, Luis Felipe Ventorim Ferrao,
//' Antonio Augusto Franco Garcia, and Matthew Stephens. 2018.
//' Harnessing Empirical Bayes and Mendelian Segregation
//' for Genotyping Autopolyploids from Messy Sequencing Data."
//' \emph{bioRxiv}. Cold Spring Harbor Laboratory. doi:10.1101/281550.
//'
//' @examples
//' ## Hardy-Weinberg population with allele-frequency of 0.75.
//' ## Moderate bias and moderate overdispersion.
//' ## See how oracle misclassification error rates change as we
//' ## increase the ploidy.
//' ploidy <- 2
//' dist <- stats::dbinom(0:ploidy, ploidy, 0.75)
//' oracle_mis(n = 100, ploidy = ploidy, seq = 0.001,
//'            bias = 0.7, od = 0.01, dist = dist)
//'
//' ploidy <- 4
//' dist <- stats::dbinom(0:ploidy, ploidy, 0.75)
//' oracle_mis(n = 100, ploidy = ploidy, seq = 0.001,
//'            bias = 0.7, od = 0.01, dist = dist)
//'
//' ploidy <- 6
//' dist <- stats::dbinom(0:ploidy, ploidy, 0.75)
//' oracle_mis(n = 100, ploidy = ploidy, seq = 0.001,
//'            bias = 0.7, od = 0.01, dist = dist)
//'
//' @export
//'
// [[Rcpp::export]]
double oracle_mis(int n,
                  int ploidy,
                  double seq,
                  double bias,
                  double od,
                  NumericVector dist) {

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


//' Returns the oracle misclassification rates for each genotype.
//'
//' Given perfect knowledge of the data generating parameters,
//' \code{oracle_mis_vec} calculates the misclassification error
//' rate at each genotype. This differs from \code{\link{oracle_mis}}
//' in that this will \emph{not} average over the genotype distribution to
//' get an overall misclassification error rate. That is, \code{oracle_mis_vec}
//' returns a vector of misclassification error rates \emph{conditional} on
//' each genotype.
//'
//' This is an ideal level of the misclassification error rate and
//' any real method will have a larger rate than this. This is a useful
//' approximation when you have a lot of individuals.
//'
//' To come up with \code{dist}, you need some additional assumptions.
//' For example, if the population is in Hardy-Weinberg equilibrium and
//' the allele frequency is \code{alpha} then you could calculate
//' \code{dist} using the R code: \code{dbinom(x = 0:ploidy, size = ploidy, prob = alpha)}.
//' Alternatively, if you know the genotypes of the individual's two parents are, say,
//' \code{ref_count1} and \code{ref_count2}, then you could use the \code{\link[updog]{get_q_array}}
//' function from the updog package: \code{get_q_array(ploidy)[ref_count1 + 1, ref_count2 + 1, ]}.
//'
//'
//' @inheritParams oracle_mis
//'
//' @return A vector of numerics. Element i is the oracle misclassification
//'     error rate when genotyping an individual with actual
//'     genotype i + 1.
//'
//' @export
//'
//' @references Gerard, David, Luis Felipe Ventorim Ferrao,
//' Antonio Augusto Franco Garcia, and Matthew Stephens. 2018.
//' Harnessing Empirical Bayes and Mendelian Segregation
//' for Genotyping Autopolyploids from Messy Sequencing Data."
//' \emph{bioRxiv}. Cold Spring Harbor Laboratory. doi:10.1101/281550.
//'
//' @examples
//' ## Hardy-Weinberg population with allele-frequency of 0.75.
//' ## Moderate bias and moderate overdispersion.
//' ploidy <- 4
//' dist <- stats::dbinom(0:ploidy, ploidy, 0.75)
//' om <- oracle_mis_vec(n = 100, ploidy = ploidy, seq = 0.001,
//'                      bias = 0.7, od = 0.01, dist = dist)
//' om
//'
//' ## Get same output as oracle_mis this way:
//' sum(dist * om)
//' oracle_mis(n = 100, ploidy = ploidy, seq = 0.001,
//'            bias = 0.7, od = 0.01, dist = dist)
//'
//' @author David Gerard
// [[Rcpp::export]]
NumericVector oracle_mis_vec(int n,
                             int ploidy,
                             double seq,
                             double bias,
                             double od,
                             NumericVector dist) {
  if (std::abs(Rcpp::sum(dist) - 1) > TOL) {
    Rcpp::stop("oracle_miss_vec: elements in dist must sum to one.");
  }
  if (dist.length() != ploidy + 1) {
    Rcpp::stop("oracle_miss_vec: dist must have the same length as ploidy + 1.");
  }

  // Calculate mean values
  NumericVector xi_vec(ploidy + 1);
  for (int i = 0; i <= ploidy; i++) {
    xi_vec(i) = xi_double((double)i / (double)ploidy, seq, bias);
  }

  // define important vectors
  // contains the misclassification errors.
  // Will be success rate until do a 1.0 - exp(mis_vec) at end.
  NumericVector mis_vec(ploidy + 1, R_NegInf);
  NumericVector ldist = Rcpp::log(dist);
  NumericVector llike(ploidy + 1); // contains the log-likelihoods

  int max_el = 0;
  double max_val = R_NegInf;
  double post_val;
  for (int x = 0; x <= n; x++) {
    // Get log-likelihoods and find maximum value
    max_el = 0;
    max_val = R_NegInf;
    for (int k = 0; k <= ploidy; k++) {
      llike(k) = dbetabinom_double(x, n, xi_vec(k), od, true);
      post_val = llike(k) + ldist(k);
      if (max_val < post_val) {
        max_val = post_val;
        max_el = k;
      }
    }

    // Add contribution to success rate.
    mis_vec(max_el) = log_sum_exp_2(mis_vec(max_el), llike(max_el));

  }

  mis_vec = 1.0 - Rcpp::exp(mis_vec);

  return mis_vec;
}


//' The joint probability of the genotype and the genotype estimate
//' of an oracle estimator.
//'
//' This returns the joint distribution of the true genotypes and an oracle
//' estimator given perfect knowledge of the data generating process. This is a useful
//' approximation when you have a lot of individuals.
//'
//' To come up with \code{dist}, you need some additional assumptions.
//' For example, if the population is in Hardy-Weinberg equilibrium and
//' the allele frequency is \code{alpha} then you could calculate
//' \code{dist} using the R code: \code{dbinom(x = 0:ploidy, size = ploidy, prob = alpha)}.
//' Alternatively, if you know the genotypes of the individual's two parents are, say,
//' \code{ref_count1} and \code{ref_count2}, then you could use the \code{\link[updog]{get_q_array}}
//' function from the updog package: \code{get_q_array(ploidy)[ref_count1 + 1, ref_count2 + 1, ]}.
//'
//' See the Examples to see how to reconcile the output of \code{oracle_joint}
//' with \code{\link{oracle_mis}} and \code{\link{oracle_mis_vec}}.
//'
//' @inheritParams oracle_mis
//'
//' @return A matrix. Element (i, j) is the joint probability of estimating
//'     the genotype to be i+1 when the true genotype is j+1. That is, the
//'     estimated genotype indexes the rows and the true genotype indexes
//'     the columns. This is when
//'     using an oracle estimator.
//'
//' @references Gerard, David, Luis Felipe Ventorim Ferrao,
//' Antonio Augusto Franco Garcia, and Matthew Stephens. 2018.
//' Harnessing Empirical Bayes and Mendelian Segregation
//' for Genotyping Autopolyploids from Messy Sequencing Data."
//' \emph{bioRxiv}. Cold Spring Harbor Laboratory. doi:10.1101/281550.
//'
//' @examples
//' ## Hardy-Weinberg population with allele-frequency of 0.75.
//' ## Moderate bias and moderate overdispersion.
//' ploidy <- 4
//' dist <- stats::dbinom(0:ploidy, ploidy, 0.75)
//' jd <- oracle_joint(n = 100, ploidy = ploidy, seq = 0.001,
//'                    bias = 0.7, od = 0.01, dist = dist)
//' jd
//'
//' ## Get same output as oracle_mis this way:
//' 1 - sum(diag(jd))
//' oracle_mis(n = 100, ploidy = ploidy, seq = 0.001,
//'            bias = 0.7, od = 0.01, dist = dist)
//'
//' ## Get same output as oracle_mis_vec this way:
//' 1 - diag(sweep(x = jd, MARGIN = 2, STATS = colSums(jd), FUN = "/"))
//' oracle_mis_vec(n = 100, ploidy = ploidy, seq = 0.001,
//'                bias = 0.7, od = 0.01, dist = dist)
//'
//' @export
//'
//' @seealso
//' \describe{
//'   \item{\code{\link{oracle_plot}}}{For visualizing the joint distribution output from \code{oracle_joint}.}
//'   \item{\code{\link{oracle_mis_from_joint}}}{For obtaining the same results as \code{\link{oracle_mis}}
//'       directly from the output of \code{oracle_joint}.}
//'   \item{\code{\link{oracle_mis_vec_from_joint}}}{For obtaining the same results as \code{\link{oracle_mis_vec}}
//'       directly from the output of \code{oracle_joint}.}
//'   \item{\code{\link{oracle_cor_from_joint}}}{For obtaining the same results as \code{\link{oracle_cor}}
//'       directly from the output of \code{oracle_joint}.}
//' }
//'
//' @author David Gerard
// [[Rcpp::export]]
NumericMatrix oracle_joint(int n,
                           int ploidy,
                           double seq,
                           double bias,
                           double od,
                           NumericVector dist) {
  if (std::abs(Rcpp::sum(dist) - 1) > TOL) {
    Rcpp::stop("oracle_miss_vec: elements in dist must sum to one.");
  }
  if (dist.length() != ploidy + 1) {
    Rcpp::stop("oracle_miss_vec: dist must have the same length as ploidy + 1.");
  }

  // Calculate mean values
  NumericVector xi_vec(ploidy + 1);
  for (int i = 0; i <= ploidy; i++) {
    xi_vec(i) = xi_double((double)i / (double)ploidy, seq, bias);
  }

  // define important vectors
  NumericMatrix joint_dist(ploidy + 1, ploidy + 1); // contains the joint distribution
  std::fill(joint_dist.begin(), joint_dist.end(), R_NegInf);
  NumericVector ldist = Rcpp::log(dist);
  NumericVector llike(ploidy + 1); // contains the log-likelihoods

  int max_el = 0;
  double max_val = R_NegInf;
  double post_val;
  for (int x = 0; x <= n; x++) {
    // Get log-likelihoods and find maximum value
    max_el = 0;
    max_val = R_NegInf;
    for (int k = 0; k <= ploidy; k++) {
      llike(k) = dbetabinom_double(x, n, xi_vec(k), od, true);
      post_val = llike(k) + ldist(k);
      if (max_val < post_val) {
        max_val = post_val;
        max_el = k;
      }
    }

    // Add contribution to joint probability
    for (int k = 0; k <= ploidy; k++) {
      joint_dist(max_el, k) = log_sum_exp_2(joint_dist(max_el, k), llike(k) + ldist(k));
    }
  }

  for (int k1 = 0; k1 <= ploidy; k1++) {
    for (int k2 = 0; k2 <= ploidy; k2++) {
      joint_dist(k1, k2) = std::exp(joint_dist(k1, k2));
    }
  }

  return joint_dist;
}
