#include "mupdog.h"

// Code to calculate the penalized ELBO.

//' Variational posterior probability of having \code{dosage} A alleles
//' when the ploidy is \code{ploidy}, the allele frequency is
//' \code{alpha}, the individual-specific overdispersion parameter is
//' \code{rho}, the variational mean is \code{mu}, and the variational
//' variance is \code{sigma2}.
//'
//' @param dosage The number of A alleles.
//' @param ploidy The ploidy of the individual.
//' @param mu The variational mean.
//' @param sigma2 The variational variance (not standard devation).
//' @param alpha The allele frequency.
//' @param rho The individual's overdispersion parameter.
//'
//' @author David
//'
// [[Rcpp::export]]
double post_prob(int dosage, int ploidy, double mu, double sigma2,
                 double alpha, double rho) {

  double pk   = pbetabinom_double(dosage, ploidy, alpha, rho, false);
  double pkm1 = pbetabinom_double(dosage - 1, ploidy, alpha, rho, false);

  double lhs = R::pnorm5(R::qnorm5(pk, 0, 1, true, false), mu,
       std::sqrt(sigma2), true, false);
  double rhs = R::pnorm5(R::qnorm5(pkm1, 0, 1, true, false), mu,
       std::sqrt(sigma2), true, false);

  double post_prob = lhs - rhs;
  return post_prob;
}

//' Computes every posterior probability for each dosage level for
//' each individual at each SNP.
//'
//' @param ploidy The ploidy of the species.
//' @param mu A matrix of variational posterior means. The rows
//'     index the individuals and the columns index the SNPs.
//' @param sigma2 A matrix of variational posterior variances.
//'     The rows index the individuals and the columns index the SNPs.
//' @param alpha A vector of allele frequencies for all SNPs.
//' @param rho A vector of inbreeding coefficients for all individuals.
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
arma::Cube<double> compute_all_post_prob(int ploidy,
                                         NumericMatrix mu,
                                         NumericMatrix sigma2,
                                         NumericVector alpha,
                                         NumericVector rho) {
  // check input ------------------------------------
  int nind = mu.nrow();
  int nsnps = mu.ncol();

  if (sigma2.nrow() != nind) {
    Rcpp::stop("sigma2 and mu must have the same number of rows.");
  }
  if (sigma2.ncol() != nsnps) {
    Rcpp::stop("sigma2 and mu must have the same number of columns.");
  }
  if (alpha.length() != nsnps) {
    Rcpp::stop("alpha must have the same length as the number of columns in mu.");
  }
  if (rho.length() != nind) {
    Rcpp::stop("rho must have the same length as the number of rows in mu.");
  }

  // iterate through all i,j,k ----------------------
  arma::Cube<double> warray(nind, nsnps, ploidy + 1);
  for (int i = 0; i < nind; i++) {
    for (int j = 0; j < nsnps; j++) {
      for (int k = 0; k <= ploidy; k++) {
        warray(i, j, k) = post_prob(k, ploidy,
               mu(i, j), sigma2(i, j), alpha(j), rho(i));
      }
    }
  }
  return warray;
}


//' Penalty on bias parameter.
//'
//' @param h The current value of bias parameter. Must be
//'     greater than 0. A value of 1 means no bias.
//' @param mu_h The prior mean of the log-bias parameter.
//' @param sigma2_h The prior variance (not standard deviation)
//'     of the log-bias parameter.
//'
//' @author David Gerard
// [[Rcpp::export]]
double pen_bias(double h, double mu_h, double sigma2_h) {

  if (h < TOL) {
    Rcpp::stop("h must be greater than 0.");
  }
  if (sigma2_h < TOL) {
    Rcpp::stop("sigma2_h must be greater tha 0.");
  }

  double pen = -std::log(h) - std::pow(std::log(h) - mu_h, 2) / (2.0 * sigma2_h);

  return pen;
}

//' Penalty on sequencing error rate.
//'
//' @param eps The current value of sequencing error rate.
//'     Must be between 0 and 1.
//' @param mu_eps The prior mean of the logit sequencing error rate.
//' @param sigma2_eps The prior variance (not standard deviation)
//'     of the logit sequencing error rate.
//'
//' @author David Gerard
// [[Rcpp::export]]
double pen_seq_error(double eps, double mu_eps, double sigma2_eps) {

  if ((eps < TOL) | ((1.0 - eps) < TOL)) {
    Rcpp::stop("eps must be between 0 and 1.");
  }
  if (sigma2_eps < TOL) {
    Rcpp::stop("sigma2_eps must be greater tha 0.");
  }
  double pen = -log(eps * (1.0 - eps)) - std::pow(logit(eps) - mu_eps, 2) / (2.0 * sigma2_eps);
  return pen;
}


//' Objective function when updating a single inbreeding coefficient.
//'
//' @param mu A vector of posterior means. The jth element is the
//'     posterior mean of SNP j for the individual.
//' @param sigma2 A vector of posterior variances. The jth element
//'     is the posterior variance of SNP j for the individual.
//' @param alpha A vector of allele frequencies. The jth element
//'     is the allele frequency for SNP j.
//' @param log_bb_dense A matrix of log posterior densities. The
//'     rows index the dosage and the columns index the SNPs.
//' @param ploidy The ploidy of the species.
//'
//'
//' @author David Gerard
// [[Rcpp::export]]
double obj_for_rho(NumericVector mu,
                   NumericVector sigma2,
                   NumericVector alpha,
                   NumericMatrix log_bb_dense,
                   int ploidy) {
  // check input ---------------------------------------


  return 0;
}







