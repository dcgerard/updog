#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

extern double TOL;

// from util_funs.cpp -----------------------------------------------
double log_sum_exp(Rcpp::NumericVector x);
double logit(double x);
double xi_double(double p, double eps, double h);

// from betabinom.cpp -----------------------------------------------
double pbetabinom_double(int q, int size, double mu, double rho, bool log_p);
double dbetabinom_double(int x, int size, double mu, double rho, bool log);
