// [[Rcpp::interfaces(r, cpp)]]

#include "mupdog.h"

double TOL = 10.0 * DBL_EPSILON; // defined here, used everywhere

//' Density function of betabinomial with the shape parameterizations
//'
//' @inheritParams dbetabinom_double
//' @param alpha The first shape parameter.
//' @param beta The second shape parameter.
//'
//' @return The density of the beta-binomial.
//'
//' @keywords internal
//' @noRd
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
//' @return The density of the Bernoulli-binomial.
//'
//' @keywords internal
//' @noRd
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
//' @return The density of the beta-binomial.
//'
//' @keywords internal
//' @noRd
//'
//' @author David Gerard
// [[Rcpp::export]]
double dbetabinom_double(int x, int size, double mu, double rho, bool log) {

  // check input --------------------------------------------
  if (size < 0) {
    Rcpp::Rcout << size << std::endl;
    Rcpp::stop("dbetabinom_double: size must be greater than 0.");
  }
  if ((x < 0) | (x > size)) {
    Rcpp::Rcout << x << std::endl;
    Rcpp::stop("dbetabinom_double: x must be between 0 and size.");
  }
  if ((mu < -TOL) | ((1.0 - mu) < -TOL)) {
    Rcpp::Rcout << mu << std::endl;
    Rcpp::stop("dbetabinom_double: mu must be between 0 and 1.");
  }
  if ((rho < -TOL) | ((1.0 - rho) < -TOL)) {
    Rcpp::Rcout << rho << std::endl;
    Rcpp::stop("dbetabinom_double: rho must be between 0 and 1.");
  }

  // calculate density.
  double dout;
  if ((mu < TOL) | ((1.0 - mu) < TOL)) {
    dout = dbernbinom(x, size, mu, log);
  }
  else if (rho < TOL) {
    dout = R::dbinom(x, size, mu, log);
  }
  else if (TOL < (1.0 - rho)) {
    double alpha = mu * (1.0 - rho) / rho;
    double beta  = (1.0 - mu) * (1.0 - rho) / rho;
    dout = dbetabinom_alpha_beta_double(x, size, alpha, beta, log);
  }
  else {
    dout = dbernbinom(x, size, mu, log);
  }
  return dout;
}


// documentation in mupdog.R
//' @describeIn betabinom Density function.
//'
//' @export
//'
// [[Rcpp::export]]
NumericVector dbetabinom(IntegerVector x, IntegerVector size,
                         NumericVector mu, NumericVector rho,
                         LogicalVector log) {
  // Check input ------------------------------------------

  int n = x.length();

  if ((n != size.length()) && (1 != size.length())) {
    Rcpp::stop("dbetabinom: size must either be of length 1 or the same length as x.");
  }
  if ((n != mu.length()) && (1 != mu.length())) {
    Rcpp::stop("dbetabinom: mu must either be of length 1 or the same length as x.");
  }
  if ((n != rho.length()) && (1 != rho.length())) {
    Rcpp::stop("dbetabinom: rho must either be of length 1 or the same length as x.");
  }
  if ((n != log.length()) && (1 != log.length())) {
    Rcpp::stop("dbetabinom: log must either be of length 1 or the same length as x.");
  }

  // iterate
  NumericVector dout(n);
  int current_size;
  double current_mu;
  double current_rho;
  bool current_log;
  for (int i = 0; i < n; i++) {
    if (size.length() == 1) {
      current_size = size(0);
    }
    else {
      current_size = size(i);
    }

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

    dout(i) = dbetabinom_double(x(i), current_size, current_mu,
         current_rho, current_log);
  }

  return dout;
}

//' The distribution function of the betabinomial. This is generally
//' only advisable if q is relatively small.
//'
//' @inheritParams dbetabinom_double
//' @param q A quantile.
//' @param log_p A logical. Should return the log-probability
//'     \code{TRUE} or not \code{FALSE}?
//'
//' @return The tail-probability of the beta-binomial.
//'
//' @keywords internal
//' @noRd
//'
//' @author David Gerard
// [[Rcpp::export]]
double pbetabinom_double(int q, int size, double mu, double rho, bool log_p) {

  // Check input ------------------------------------------------------
  if (size < 0) {
    Rcpp::Rcout << size << std::endl;
    Rcpp::stop("pbetabinom_double: size must be greater than 0.");
  }
  if ((mu < -TOL) | ((1.0 - mu) < -TOL)) {
    Rcpp::Rcout << mu << std::endl;
    Rcpp::stop("pbetabinom_double: mu must be between 0 and 1.");
  }
  if ((rho < -TOL) | ((1.0 - rho) < -TOL)) {
    Rcpp::Rcout << rho << std::endl;
    Rcpp::stop("pbetabinom_double: rho must be between 0 and 1.");
  }

  double lp; // the log of the cdf.
  if (q > size) {
    lp = 0.0;
  }
  else if (q < 0) {
    lp = R_NegInf;
  }
  else {
    NumericVector log_prob_vec(q + 1);
    for (int i = 0; i <= q; i++) {
      log_prob_vec(i) = dbetabinom_double(i, size, mu, rho, true);
    }
    lp = log_sum_exp(log_prob_vec);
  }

  if (log_p) {
    return lp;
  }
  else {
    return std::min(std::exp(lp), 1.0);
  }
}

// documentation in mupdog.R
//' @describeIn betabinom Distribution function.
//'
//' @export
//'
// [[Rcpp::export]]
NumericVector pbetabinom(IntegerVector q, IntegerVector size,
                         NumericVector mu, NumericVector rho,
                         LogicalVector log_p) {

  // Check input ------------------------------------------

  int n = q.length();

  if ((n != size.length()) && (1 != size.length())) {
    Rcpp::stop("pbetabinom: size must either be of length 1 or the same length as q.");
  }
  if ((n != mu.length()) && (1 != mu.length())) {
    Rcpp::stop("pbetabinom: mu must either be of length 1 or the same length as q.");
  }
  if ((n != rho.length()) && (1 != rho.length())) {
    Rcpp::stop("pbetabinom: rho must either be of length 1 or the same length as q.");
  }
  if ((n != log_p.length()) && (1 != log_p.length())) {
    Rcpp::stop("pbetabinom: log_p must either be of length 1 or the same length as q.");
  }

  // iterate
  NumericVector pout(n);
  int current_size;
  double current_mu;
  double current_rho;
  bool current_log_p;
  for (int i = 0; i < n; i++) {
    if (size.length() == 1) {
      current_size = size(0);
    }
    else {
      current_size = size(i);
    }

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

    if (log_p.length() == 1) {
      current_log_p = log_p(0);
    }
    else {
      current_log_p = log_p(i);
    }

    pout(i) = pbetabinom_double(q(i), current_size, current_mu,
         current_rho, current_log_p);
  }

  return pout;
}


//' The quantile function of the beta-binomial distribution parameterized
//' by mean and overdispersion parameter.
//'
//' @inheritParams dbetabinom_double
//' @param p The lower tail probability.
//'
//' @return The quantile of the beta-binomial.
//'
//' @keywords internal
//' @noRd
//'
//' @author David Gerard
// [[Rcpp::export]]
int qbetabinom_double(double p, int size, double mu, double rho) {

  // check input -------------------------------------------------
  if ((p < 0.0) | (p > 1.0)) {
    Rcpp::stop("qbetabinom_double: p needs to be between 0 and 1.");
  }
  if ((mu < 0.0) | (mu > 1.0)) {
    Rcpp::stop("qbetabinom_double: mu needs to be between 0 and 1.");
  }
  if ((rho < 0.0) | (rho > 1.0)) {
    Rcpp::stop("qbetabinom_double: rho needs to be between 0 and 1.");
  }
  if (size < 0) {
    Rcpp::stop("qbetabinom_double: size needs to be greater than or equal to 0.");
  }


  double p_current = 0.0;
  int qfinal = -1;

  if (p > 1.0 - TOL) { // this deals with a couple numerical issues I had during unit tests.
    qfinal = size;
  }
  else {
    for (int i = 0; i <= size; i++) {
      p_current = p_current + dbetabinom_double(i, size, mu, rho, false);
      if (p_current > p - TOL) {
        qfinal = i;
        break;
      }
    }
  }


  return qfinal;
}


// documentation in updog.R
//' @describeIn betabinom Quantile function.
//'
//' @export
//'
// [[Rcpp::export]]
IntegerVector qbetabinom(NumericVector p, IntegerVector size,
                         NumericVector mu, NumericVector rho) {
  // Check input ------------------------------------------
  int n = p.length();

  if ((n != size.length()) && (1 != size.length())){
    Rcpp::stop("qbetabinom: size must be either of lenght 1 or the same length as p.");
  }
  if ((n != mu.length()) && (1 != mu.length())) {
    Rcpp::stop("qbetabinom: mu must either be of length 1 or the same length as p.");
  }
  if ((n != rho.length()) && (1 != rho.length())) {
    Rcpp::stop("qbetabinom: rho must either be of length 1 or the same length as p.");
  }

  // Fit the qbetabinoms
  IntegerVector qout(n);
  int current_size;
  double current_mu;
  double current_rho;
  for (int i = 0; i < n; i++) {

    if (size.length() == 1) {
      current_size = size(0);
    }
    else {
      current_size = size(i);
    }

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

    qout(i) = qbetabinom_double(p(i), current_size, current_mu, current_rho);
  }

  return qout;
}



//' One draw from the beta-binomial distribution parameterized
//' by mean and overdispersion parameter.
//'
//' @inheritParams dbetabinom_double
//'
//' @return A random sample from the beta-binomial.
//'
//' @keywords internal
//' @noRd
//'
//' @author David Gerard
// [[Rcpp::export]]
int rbetabinom_int(int size, double mu, double rho) {
  double alpha = mu * (1.0 - rho) / rho;
  double beta = (1.0 - mu) * (1.0 - rho) / rho;
  double currentp = R::rbeta(alpha, beta);
  int draw = R::rbinom(size, currentp);
  return draw;
}

// documentation in updog.R
//' @describeIn betabinom Random generation.
//'
//' @export
//'
// [[Rcpp::export]]
IntegerVector rbetabinom(int n, IntegerVector size,
                         NumericVector mu, NumericVector rho) {
  // Check input ------------------------------------------
  if ((n != size.length()) && (1 != size.length())){
    Rcpp::stop("qbetabinom: size must be either of length 1 or n.");
  }
  if ((n != mu.length()) && (1 != mu.length())) {
    Rcpp::stop("qbetabinom: mu must either be of length 1 or n.");
  }
  if ((n != rho.length()) && (1 != rho.length())) {
    Rcpp::stop("qbetabinom: rho must either be of length 1 or n.");
  }

  IntegerVector draws(n);
  int current_size;
  double current_mu;
  double current_rho;
  for (int i = 0; i < n; i++) {
    if (size.length() == 1) {
      current_size = size(0);
    }
    else {
      current_size = size(i);
    }

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

    draws(i) = rbetabinom_int(current_size, current_mu, current_rho);
  }

  return draws;
}









