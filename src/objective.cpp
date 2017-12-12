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
