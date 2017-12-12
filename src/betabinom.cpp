#include "mupdog.h"

double TOL = 2.0 * DBL_EPSILON; // defined here, used everywhere

//' Density function of betabinomial with the shape parameterizations
//'
//' @inheritParams dbetabinom_double
//' @param alpha The first shape parameter.
//' @param beta The second shape paramter.
//'
//' @author David Gerard
// [[Rcpp::export]]
double dbetabinom_alpha_beta_double(int x, int size, double alpha, double beta, bool log) {
  double ldense = R::lchoose(size, x) +
    R::lbeta((double)x + alpha, (double)size - (double)x + beta) -
    R::lbeta(alpha, beta);

  if (log) {
    return(ldense);
  }
  else {
    return(std::exp(ldense));
  }
}

//' Special case of betabinomial where the beta is bernoulli mu.
//'
//' @inheritParams dbetabinom_double
//'
//' @author David Gerard
// [[Rcpp::export]]
double dbernbinom(int x, int size, double mu, bool log) {
  double dout; // the output density
  if (x == size) {
    if (mu > TOL) {
      dout = std::log(mu);
    }
    else {
      dout = R_NegInf;
    }
  }
  else if (x == 0) {
    if ((1.0 - mu) > TOL) {
      dout = std::log(1.0 - mu);
    }
    else {
      dout = R_NegInf;
    }
  }
  else {
    dout = R_NegInf;
  }

  if (log) {
    return dout;
  }
  else {
    return std::exp(dout);
  }

}

//' The density function of the beta-binomial distribution.
//'
//' @param x The quantile.
//' @param size The total number of draws.
//' @param mu The mean of the beta.
//' @param rho The overdispersion parameter of the beta.
//' @param log A logical. Should we return the log of the
//'     density \code{TRUE} or not \code{FALSE}?
//'
//' @author David Gerard
// [[Rcpp::export]]
double dbetabinom_double(int x, int size, double mu, double rho, bool log) {

  // check input --------------------------------------------
  if (size < 0) {
    Rcpp::stop("size must be greater than 0.");
  }
  if ((x < 0) | (x > size)) {
    Rcpp::stop("x must be between 0 and size.");
  }
  if ((mu < -TOL) | ((1.0 - mu) < -TOL)) {
    Rcpp::stop("mu must be between 0 and 1.");
  }
  if ((rho < -TOL) | ((1.0 - rho) < -TOL)) {
    Rcpp::stop("rho must be between 0 and 1.");
  }

  // calculate density.
  double dout;
  if ((mu < TOL) | ((1.0 - mu) < TOL)) {
    dout = dbernbinom(x, size, mu, log);
  }
  else if (rho < TOL) {
    dout = R::dbinom(x, size, mu, log);
  }
  else if (TOL < 1.0 - rho) {
    double alpha = mu * (1.0 - rho) / rho;
    double beta  = (1.0 - mu) * (1.0 - rho) / rho;
    dout = dbetabinom_alpha_beta_double(x, size, alpha, beta, log);
  }
  else {
    dout = dbernbinom(x, size, mu, log);
  }
  return dout;
}

//' Density function for the betabinomial distribution parameterized
//' by the mean and the overdispersion parameter.
//'
//' This is the density function. It handles the corner cases when
//' rho is 0 or 1 in intelligible ways.
//'
//' Let \eqn{\mu} and \eqn{\rho} be the mean and overdispersion paramters.
//' Let \eqn{\alpha} and \eqn{\beta} be the usual shape parameters of
//' a beta distribution. Then we have the relation
//' \deqn{\mu = \alpha/(\alpha + \beta),}
//' and
//' \deqn{\rho = 1/(1 + \alpha + \beta).}
//' This necessarily means that
//' \deqn{\alpha = \mu (1 - \rho)/\rho,}
//' and
//' \deqn{\beta = (1 - \mu) (1 - \rho)/\rho.}
//'
//' @param x A vector of quantiles.
//' @param size A vector of sizes.
//' @param mu Either a scalar of the mean for each observation,
//'     or a vector of means of each observation, and thus
//'     the same length as \code{x} and \code{size}. This must
//'     be between 0 and 1.
//' @param rho Either a scalar of the overdispersion parameter
//'     for each observation, or a vector of overdispersion
//'     parameters of each observation, and thus the same length as
//'     \code{x} and \code{size}. This must be between 0 and 1.
//' @param log A logical vector either of length 1 or the same
//'     length as \code{x} and \code{size}. This determines whether
//'     to return the log density for all observations (in the case
//'     that its length is 1) or for each observation (in the case that
//'     its length is that of \code{x} and \code{size}).
//'
//' @author David Gerard
// [[Rcpp::export]]
NumericVector dbetabinom(IntegerVector x, IntegerVector size,
                         NumericVector mu, NumericVector rho,
                         LogicalVector log) {
  // Check input ------------------------------------------

  int n = x.length();

  if (n != size.length()) {
    Rcpp::stop("x and size must be of same length.");
  }
  if ((n != mu.length()) & (1 != mu.length())) {
    Rcpp::stop("mu must either be of length 1 or the same length as x.");
  }
  if ((n != rho.length()) & (1 != rho.length())) {
    Rcpp::stop("rho must either be of length 1 or the same length as x.");
  }
  if ((n != log.length()) & (1 != log.length())) {
    Rcpp::stop("log must either be of length 1 or the same length as x.");
  }

  // iterate
  NumericVector dout(n);
  double current_mu;
  double current_rho;
  bool current_log;
  for (int i = 0; i < n; i++) {
    if (mu.length() == 1) {
      current_mu = mu(0);
    }
    else {
      current_mu = mu(i);
    }

    if (rho.length() == 1) {
      current_rho = rho(0);
    }
    else {
      current_rho = rho(i);
    }

    if (log.length() == 1) {
      current_log = log(0);
    }
    else {
      current_log = log(i);
    }

    dout(i) = dbetabinom_double(x(i), size(i), current_mu,
         current_rho, current_log);
  }

  return dout;
}
