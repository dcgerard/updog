// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// dbetabinom_alpha_beta_double
double dbetabinom_alpha_beta_double(int x, int size, double alpha, double beta, bool log);
RcppExport SEXP _mupdog_dbetabinom_alpha_beta_double(SEXP xSEXP, SEXP sizeSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< bool >::type log(logSEXP);
    rcpp_result_gen = Rcpp::wrap(dbetabinom_alpha_beta_double(x, size, alpha, beta, log));
    return rcpp_result_gen;
END_RCPP
}
// dbernbinom
double dbernbinom(int x, int size, double mu, bool log);
RcppExport SEXP _mupdog_dbernbinom(SEXP xSEXP, SEXP sizeSEXP, SEXP muSEXP, SEXP logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< bool >::type log(logSEXP);
    rcpp_result_gen = Rcpp::wrap(dbernbinom(x, size, mu, log));
    return rcpp_result_gen;
END_RCPP
}
// dbetabinom_double
double dbetabinom_double(int x, int size, double mu, double rho, bool log);
RcppExport SEXP _mupdog_dbetabinom_double(SEXP xSEXP, SEXP sizeSEXP, SEXP muSEXP, SEXP rhoSEXP, SEXP logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< bool >::type log(logSEXP);
    rcpp_result_gen = Rcpp::wrap(dbetabinom_double(x, size, mu, rho, log));
    return rcpp_result_gen;
END_RCPP
}
// dbetabinom
NumericVector dbetabinom(IntegerVector x, IntegerVector size, NumericVector mu, NumericVector rho, LogicalVector log);
RcppExport SEXP _mupdog_dbetabinom(SEXP xSEXP, SEXP sizeSEXP, SEXP muSEXP, SEXP rhoSEXP, SEXP logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type log(logSEXP);
    rcpp_result_gen = Rcpp::wrap(dbetabinom(x, size, mu, rho, log));
    return rcpp_result_gen;
END_RCPP
}
// pbetabinom_double
double pbetabinom_double(int q, int size, double mu, double rho, bool log_p);
RcppExport SEXP _mupdog_pbetabinom_double(SEXP qSEXP, SEXP sizeSEXP, SEXP muSEXP, SEXP rhoSEXP, SEXP log_pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< bool >::type log_p(log_pSEXP);
    rcpp_result_gen = Rcpp::wrap(pbetabinom_double(q, size, mu, rho, log_p));
    return rcpp_result_gen;
END_RCPP
}
// pbetabinom
NumericVector pbetabinom(IntegerVector q, IntegerVector size, NumericVector mu, NumericVector rho, LogicalVector log_p);
RcppExport SEXP _mupdog_pbetabinom(SEXP qSEXP, SEXP sizeSEXP, SEXP muSEXP, SEXP rhoSEXP, SEXP log_pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type log_p(log_pSEXP);
    rcpp_result_gen = Rcpp::wrap(pbetabinom(q, size, mu, rho, log_p));
    return rcpp_result_gen;
END_RCPP
}
// eta_double
double eta_double(double p, double eps);
RcppExport SEXP _mupdog_eta_double(SEXP pSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(eta_double(p, eps));
    return rcpp_result_gen;
END_RCPP
}
// eta_fun
NumericVector eta_fun(NumericVector p, double eps);
RcppExport SEXP _mupdog_eta_fun(SEXP pSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(eta_fun(p, eps));
    return rcpp_result_gen;
END_RCPP
}
// xi_fun
NumericVector xi_fun(NumericVector p, double eps, double h);
RcppExport SEXP _mupdog_xi_fun(SEXP pSEXP, SEXP epsSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(xi_fun(p, eps, h));
    return rcpp_result_gen;
END_RCPP
}
// log_sum_exp
double log_sum_exp(NumericVector x);
RcppExport SEXP _mupdog_log_sum_exp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(log_sum_exp(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mupdog_dbetabinom_alpha_beta_double", (DL_FUNC) &_mupdog_dbetabinom_alpha_beta_double, 5},
    {"_mupdog_dbernbinom", (DL_FUNC) &_mupdog_dbernbinom, 4},
    {"_mupdog_dbetabinom_double", (DL_FUNC) &_mupdog_dbetabinom_double, 5},
    {"_mupdog_dbetabinom", (DL_FUNC) &_mupdog_dbetabinom, 5},
    {"_mupdog_pbetabinom_double", (DL_FUNC) &_mupdog_pbetabinom_double, 5},
    {"_mupdog_pbetabinom", (DL_FUNC) &_mupdog_pbetabinom, 5},
    {"_mupdog_eta_double", (DL_FUNC) &_mupdog_eta_double, 2},
    {"_mupdog_eta_fun", (DL_FUNC) &_mupdog_eta_fun, 2},
    {"_mupdog_xi_fun", (DL_FUNC) &_mupdog_xi_fun, 3},
    {"_mupdog_log_sum_exp", (DL_FUNC) &_mupdog_log_sum_exp, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_mupdog(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
