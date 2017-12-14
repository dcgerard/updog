// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
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
// post_prob
double post_prob(int dosage, int ploidy, double mu, double sigma2, double alpha, double rho);
RcppExport SEXP _mupdog_post_prob(SEXP dosageSEXP, SEXP ploidySEXP, SEXP muSEXP, SEXP sigma2SEXP, SEXP alphaSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type dosage(dosageSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(post_prob(dosage, ploidy, mu, sigma2, alpha, rho));
    return rcpp_result_gen;
END_RCPP
}
// compute_all_post_prob
arma::Cube<double> compute_all_post_prob(int ploidy, NumericMatrix mu, NumericMatrix sigma2, NumericVector alpha, NumericVector rho);
RcppExport SEXP _mupdog_compute_all_post_prob(SEXP ploidySEXP, SEXP muSEXP, SEXP sigma2SEXP, SEXP alphaSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_all_post_prob(ploidy, mu, sigma2, alpha, rho));
    return rcpp_result_gen;
END_RCPP
}
// compute_all_log_bb
arma::Cube<double> compute_all_log_bb(NumericMatrix refmat, NumericMatrix sizemat, int ploidy, NumericVector seq, NumericVector bias, NumericVector od);
RcppExport SEXP _mupdog_compute_all_log_bb(SEXP refmatSEXP, SEXP sizematSEXP, SEXP ploidySEXP, SEXP seqSEXP, SEXP biasSEXP, SEXP odSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type refmat(refmatSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type sizemat(sizematSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bias(biasSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type od(odSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_all_log_bb(refmat, sizemat, ploidy, seq, bias, od));
    return rcpp_result_gen;
END_RCPP
}
// pen_bias
double pen_bias(double h, double mu_h, double sigma2_h);
RcppExport SEXP _mupdog_pen_bias(SEXP hSEXP, SEXP mu_hSEXP, SEXP sigma2_hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type mu_h(mu_hSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2_h(sigma2_hSEXP);
    rcpp_result_gen = Rcpp::wrap(pen_bias(h, mu_h, sigma2_h));
    return rcpp_result_gen;
END_RCPP
}
// pen_seq_error
double pen_seq_error(double eps, double mu_eps, double sigma2_eps);
RcppExport SEXP _mupdog_pen_seq_error(SEXP epsSEXP, SEXP mu_epsSEXP, SEXP sigma2_epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< double >::type mu_eps(mu_epsSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2_eps(sigma2_epsSEXP);
    rcpp_result_gen = Rcpp::wrap(pen_seq_error(eps, mu_eps, sigma2_eps));
    return rcpp_result_gen;
END_RCPP
}
// obj_for_rho
double obj_for_rho(double rho, NumericVector mu, NumericVector sigma2, NumericVector alpha, NumericMatrix log_bb_dense, int ploidy);
RcppExport SEXP _mupdog_obj_for_rho(SEXP rhoSEXP, SEXP muSEXP, SEXP sigma2SEXP, SEXP alphaSEXP, SEXP log_bb_denseSEXP, SEXP ploidySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type log_bb_dense(log_bb_denseSEXP);
    Rcpp::traits::input_parameter< int >::type ploidy(ploidySEXP);
    rcpp_result_gen = Rcpp::wrap(obj_for_rho(rho, mu, sigma2, alpha, log_bb_dense, ploidy));
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
NumericVector eta_fun(NumericVector p, NumericVector eps);
RcppExport SEXP _mupdog_eta_fun(SEXP pSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(eta_fun(p, eps));
    return rcpp_result_gen;
END_RCPP
}
// xi_double
double xi_double(double p, double eps, double h);
RcppExport SEXP _mupdog_xi_double(SEXP pSEXP, SEXP epsSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(xi_double(p, eps, h));
    return rcpp_result_gen;
END_RCPP
}
// xi_fun
NumericVector xi_fun(NumericVector p, NumericVector eps, NumericVector h);
RcppExport SEXP _mupdog_xi_fun(SEXP pSEXP, SEXP epsSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
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
// logit
double logit(double x);
RcppExport SEXP _mupdog_logit(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(logit(x));
    return rcpp_result_gen;
END_RCPP
}
// expit
double expit(double x);
RcppExport SEXP _mupdog_expit(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(expit(x));
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
    {"_mupdog_post_prob", (DL_FUNC) &_mupdog_post_prob, 6},
    {"_mupdog_compute_all_post_prob", (DL_FUNC) &_mupdog_compute_all_post_prob, 5},
    {"_mupdog_compute_all_log_bb", (DL_FUNC) &_mupdog_compute_all_log_bb, 6},
    {"_mupdog_pen_bias", (DL_FUNC) &_mupdog_pen_bias, 3},
    {"_mupdog_pen_seq_error", (DL_FUNC) &_mupdog_pen_seq_error, 3},
    {"_mupdog_obj_for_rho", (DL_FUNC) &_mupdog_obj_for_rho, 6},
    {"_mupdog_eta_double", (DL_FUNC) &_mupdog_eta_double, 2},
    {"_mupdog_eta_fun", (DL_FUNC) &_mupdog_eta_fun, 2},
    {"_mupdog_xi_double", (DL_FUNC) &_mupdog_xi_double, 3},
    {"_mupdog_xi_fun", (DL_FUNC) &_mupdog_xi_fun, 3},
    {"_mupdog_log_sum_exp", (DL_FUNC) &_mupdog_log_sum_exp, 1},
    {"_mupdog_logit", (DL_FUNC) &_mupdog_logit, 1},
    {"_mupdog_expit", (DL_FUNC) &_mupdog_expit, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_mupdog(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
