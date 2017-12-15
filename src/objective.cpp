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

  if ((alpha < -TOL) | ((1.0 - alpha) < -TOL)) {
    Rcpp::Rcout << alpha << std::endl;
    Rcpp::stop("post_prob: alpha must be between 0 and 1.");
  }
  if ((rho < -TOL) | ((1.0 - rho) < -TOL)) {
    Rcpp::Rcout << rho << std::endl;
    Rcpp::stop("post_prob: rho must be between 0 and 1.");
  }
  if (sigma2 < -TOL) {
    Rcpp::Rcout << sigma2 << std::endl;
    Rcpp::stop("post_prob: sigma2 must be greater than 0.");
  }

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
    Rcpp::stop("compute_all_post_prob: sigma2 and mu must have the same number of rows.");
  }
  if (sigma2.ncol() != nsnps) {
    Rcpp::stop("compute_all_post_prob: sigma2 and mu must have the same number of columns.");
  }
  if (alpha.length() != nsnps) {
    Rcpp::stop("compute_all_post_prob: alpha must have the same length as the number of columns in mu.");
  }
  if (rho.length() != nind) {
    Rcpp::stop("compute_all_post_prob: rho must have the same length as the number of rows in mu.");
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

//' Calculates the log-density for every individual by snp by dosage level.
//'
//' @inheritParams mupdog
//'
//'
//'
// [[Rcpp::export]]
arma::Cube<double> compute_all_log_bb(NumericMatrix refmat, NumericMatrix sizemat,
				      int ploidy, NumericVector seq, NumericVector bias, NumericVector od) {

  // Check input ---------------------------------------
  int nind  = refmat.nrow();
  int nsnps = refmat.ncol();
  if (sizemat.nrow() != nind) {
    Rcpp::stop("compute_all_log_bb: sizemat and refmat must have the same dimensions.");
  }
  if (sizemat.ncol() != nsnps) {
    Rcpp::stop("compute_all_log_bb: sizemat and refmat must have the same dimensions.");
  }
  if (seq.length() != nsnps) {
    Rcpp::stop("compute_all_log_bb: seq must have the same length as the number of SNPs.");
  }
  if (bias.length() != nsnps) {
    Rcpp::stop("compute_all_log_bb: bias must have the same length as the number of SNPs.");
  }
  if (od.length() != nsnps) {
    Rcpp::stop("compute_all_log_bb: od must have the same length as the number of SNPs.");
  }

  arma::Cube<double> logbbdense(nind, nsnps, ploidy + 1);

  double xi_current;
  for (int i = 0; i < nind; i++) {
    for (int j = 0; j < nsnps; j++) {
      for (int k = 0; k <= ploidy; k++) {
	if (R_IsNA(refmat(i, j)) | R_IsNA(sizemat(i, j))) {
	  logbbdense(i, j, k) = NA_REAL;
	}
	else {
	  xi_current = xi_double((double)k / (double)ploidy, seq(j), bias(j));
	  logbbdense(i, j, k) = dbetabinom_double(refmat(i, j), sizemat(i, j), xi_current, od(j), true);
	}
      }
    }
  }

  return logbbdense;
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
    Rcpp::stop("pen_bias: h must be greater than 0.");
  }
  if (sigma2_h < TOL) {
    Rcpp::stop("pen_bias: sigma2_h must be greater tha 0.");
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
    Rcpp::stop("pen_seq_error: eps must be between 0 and 1.");
  }
  if (sigma2_eps < TOL) {
    Rcpp::stop("pen_seq_error: sigma2_eps must be greater tha 0.");
  }
  double pen = -log(eps * (1.0 - eps)) - std::pow(logit(eps) - mu_eps, 2) / (2.0 * sigma2_eps);
  return pen;
}


//' Objective function when updating a single inbreeding coefficient.
//'
//' @param rho The inbreeding coefficient for the individual.
//' @param mu A vector of posterior means. The jth element is the
//'     posterior mean of SNP j for the individual.
//' @param sigma2 A vector of posterior variances. The jth element
//'     is the posterior variance of SNP j for the individual.
//' @param alpha A vector of allele frequencies. The jth element
//'     is the allele frequency for SNP j.
//' @param log_bb_dense A matrix of log posterior densities. The
//'     rows index the SNPs and the columns index the dosage.
//' @param ploidy The ploidy of the species.
//'
//'
//' @author David Gerard
// [[Rcpp::export]]
double obj_for_rho(double rho,
		   NumericVector mu,
                   NumericVector sigma2,
                   NumericVector alpha,
                   NumericMatrix log_bb_dense,
                   int ploidy) {
  // check input ---------------------------------------
  int nsnps = log_bb_dense.nrow();
  if (log_bb_dense.ncol() != ploidy + 1) {
    Rcpp::stop("obj_for_rho: log_bb_dense must have ploidy+1 columns.");
  }
  if (mu.length() != nsnps) {
    Rcpp::stop("obj_for_rho: mu must have length equal to the number of SNPs");
  }
  if (sigma2.length() != nsnps) {
    Rcpp::stop("obj_for_rho: sigma2 must have length equal to the number of SNPs");
  }
  if (alpha.length() != nsnps) {
    Rcpp::stop("obj_for_rho: alpha must have length equal to the number of SNPs");
  }

  // calculate objective -------------------------------
  double obj_val = 0.0;

  double w_current;
  for (int j = 0; j < nsnps; j++) {
    for (int k = 0; k <= ploidy; k++) {
      if (!R_IsNA(log_bb_dense(j, k))) {
	w_current = post_prob(k, ploidy, mu(j), sigma2(j),
			      alpha(j), rho);
	obj_val = obj_val + w_current * log_bb_dense(j, k);
      }
    }
  }

  return obj_val;
}


//' Objective function when updating mu, sigma2, and alpha
//'
//' @param mu A vector. The ith element is individual i's variational posterior mean at the SNP.
//' @param sigma2 A vector. The ith element is individual i's variational posterior variance at the SNP.
//' @param alpha The SNP's allele frequency.
//' @param rho A vector. The ith element is individuals i's inbreeding coefficient.
//' @param log_bb_dense A matrix of log-densities of the beta binomial. The rows index the individuals and the columns index the allele dosage.
//' @param ploidy The ploidy of the species.
//' @param cor_inv The inverse of the correlation matrix.
//'
//'
//'
//' @author David Gerard
// [[Rcpp::export]]
double obj_for_mu(arma::Col<double> mu,
		  arma::Col<double> sigma2,
		  double alpha,
		  NumericVector rho,
		  NumericMatrix log_bb_dense,
		  int ploidy,
		  arma::Mat<double> cor_inv) {

  // check input --------------------------------------------------
  int nind = log_bb_dense.nrow();
  if (log_bb_dense.ncol() != ploidy + 1) {
    Rcpp::stop("obj_for_mu: log_bb_dense must have ploidy+1 columns.");
  }
  if (mu.n_elem != nind) {
    Rcpp::stop("obj_for_mu: mu must have length equal to the number of individuals.");
  }
  if (sigma2.n_elem != nind) {
    Rcpp::stop("obj_for_mu: sigma2 must have length equal to the number of individuals.");
  }
  if (rho.length() != nind) {
    Rcpp::stop("obj_for_mu: rho must have length equal to the number of individuals.");
  }
  if (cor_inv.n_rows != nind) {
    Rcpp::stop("obj_for_mu: cor_inv must have nrow equal to the number of individuals.");
  }
  if (cor_inv.n_cols != nind) {
    Rcpp::stop("obj_for_mu: cor_inv must have ncol equal to the number of individuals.");
  }

  // likelihood integral ----------------------------------------------------
  double obj_val = 0.0;
  double w_current;
  for (int i = 0; i < nind; i++) {
    for (int k = 0; k <= ploidy; k++) {
       if (!R_IsNA(log_bb_dense(i, k))) {
	 w_current = post_prob(k, ploidy, mu(i), sigma2(i), alpha, rho(i));
	 obj_val = obj_val + w_current * log_bb_dense(i, k);
       }
    }
  }

  // Multivariate normal prior integral -------------------------------------
  arma::Mat<double> muRmu = mu.t() * cor_inv * mu;
  arma::Mat<double> trRsigma2 = arma::diagvec(cor_inv).t() * sigma2;
  obj_val = obj_val - muRmu(0, 0) / 2.0 - trRsigma2(0, 0) / 2.0 + arma::accu(arma::log(sigma2)) / 2.0;

  return obj_val;
}


//' Wrapper for \code{\link{obj_for_mu}}.
//'
//' @inheritParams obj_for_mu
//' @param muSigma2Alpha A vector where the first \code{nsnps} observations are mu,
//'     the next \code{nsnps} observations are sigma2, and the last observation is \code{alpha}.
//'
//'
//'
//'
//'
//'
//' @author David Gerard
// [[Rcpp::export]]
double obj_for_mu_wrapper(arma::Col<double> muSigma2Alpha,
			  NumericVector rho,
			  NumericMatrix log_bb_dense,
			  int ploidy,
			  arma::Mat<double> cor_inv) {
  int nsnps = log_bb_dense.nrow();

  if (muSigma2Alpha.n_elem != nsnps * 2 + 1) {
    Rcpp::stop("obj_for_mu_wrapper: muSigma2Alpha must have nrow(log_bb_dense) * 2 + 1 elements");
  }

  arma::Col<double> mu(nsnps);
  arma::Col<double> sigma2(nsnps);
  double alpha;

  for (int i = 0; i < nsnps; i++) {
    mu(i) = muSigma2Alpha(i);
  }
  for (int i = nsnps; i < (2 * nsnps); i++) {
    sigma2(i - nsnps) = muSigma2Alpha(nsnps);
  }
  alpha = muSigma2Alpha(2 * nsnps);

  // Rcpp::Rcout << alpha << std::endl;
  // if (alpha < 0) {
  //   Rcpp::Rcout << mu.t() << std::endl << sigma2.t() << std::endl;
  // }


  double obj_val = obj_for_mu(mu, sigma2, alpha, rho, log_bb_dense, ploidy, cor_inv);

  // Rcpp::Rcout << obj_val << std::endl;

  return obj_val;
}


//' Objective function for updating sequencing error rate, bias, and overdispersion parameters.
//'
//' @param parvec A vector of length three. The first element is the sequencing
//'     error rate, the second element is the allele bias, and the third element
//'     is the overdispersion parameter.
//' @param refvec A vector. The ith element is the reference count for the ith individual in the SNP.
//' @param sizevec A vector. the ith element is the size count for the ith individual in the SNP/
//' @param ploidy The ploidy of the species.
//' @param mean_bias The prior mean of the log-bias.
//' @param var_bias The prior variance of the log-bias
//' @param mean_seq The prior mean of the logit sequencing error rate.
//' @param var_seq The prior variance of the logit sequencing error rate.
//' @param wmat The matrix of variational posterior probabilities for each dosage.
//'     The rows index the individuals and the columns index the dosage levels.
//'
//' @author David Gerard
// [[Rcpp::export]]
double obj_for_eps(NumericVector parvec,
		   NumericVector refvec,
		   NumericVector sizevec,
		   int ploidy,
		   double mean_bias,
		   double var_bias,
		   double mean_seq,
		   double var_seq,
		   NumericMatrix wmat) {
  // Check input -------------------------------------------------------
  int nind = wmat.nrow();

  if (parvec.length() != 3) {
    Rcpp::Rcout << parvec.length();
    Rcpp::stop("obj_for_eps: parvec must have length 3.");
  }
  if (refvec.length() != nind) {
    Rcpp::Rcout << refvec.length();
    Rcpp::stop("obj_for_eps: refvec needs to have the same length as the number of individuals.");
  }
  if (sizevec.length() != nind) {
    Rcpp::Rcout << sizevec.length();
    Rcpp::stop("obj_for_eps: sizevec needs to have the same length as the number of individuals.");
  }
  if (wmat.ncol() != ploidy + 1) {
    Rcpp::Rcout << wmat.ncol();
    Rcpp::stop("obj_for_eps: wmat needs to have ploidy+1 columns.");
  }

  double eps = parvec(0);
  double h   = parvec(1);
  double tau = parvec(2);
  double xi_current;
  double obj_val = 0.0;

  for (int i = 0; i < nind; i++) {
    for (int k = 0; k <= ploidy; k++) {
      if (!R_IsNA(refvec(i)) & !R_IsNA(sizevec(i))) {
	xi_current = xi_double((double)k / (double)ploidy, eps, h);
	obj_val = obj_val + wmat(i, k) * dbetabinom_double(refvec(i), sizevec(i), xi_current, tau, true);
      }
    }
  }

  obj_val = obj_val + pen_bias(h, mean_bias, var_bias) + pen_seq_error(eps, mean_seq, var_seq);

  if (obj_val == R_NegInf) {
    Rcpp::Rcout << obj_val << std::endl;
    Rcpp::Rcout << parvec << std::endl;
  }

  return obj_val;
}


//' The evidence lower bound
//'
//'
//'
//'
//'
//'
//'
//'
//'
//'
//'
//'
//'
//'
//'
// [[Rcpp::export]]
double elbo(arma::Cube<double> warray, arma::Cube<double> lbeta_array,
	    arma::Mat<double> cor_inv, arma::Mat<double> postmean,
	    arma::Mat<double> postvar,
	    NumericVector bias, NumericVector seq, double mean_bias,
	    double var_bias, double mean_seq, double var_seq, int ploidy) {
  // Check input -------------------------------------------------------
  int nsnps = warray.n_cols;
  int nind  = warray.n_rows;

  double obj_val = 0.0;

  // From likelihood
  for (int i = 0; i < nind; i++) {
    for (int j = 0; j < nsnps; j++) {
      for (int k = 0; k <= ploidy; k++) {
	if (!R_IsNA(lbeta_array(i, j, k))) {
	  obj_val = obj_val + warray(i, j, k) * lbeta_array(i, j, k);
	}
      }
    }
  }

  // From prior
  obj_val = obj_val + arma::log_det(cor_inv).real() * (double)nsnps / 2.0; // add rather than subtract because using inverse correlation matrix.
  obj_val = obj_val + arma::accu(arma::log(postvar)) / 2.0;

  arma::Col<double> diag_cor_inv = diagvec(cor_inv);

  double mu_cont;
  double var_cont;
  for (int j = 0; j < nsnps; j++) {
    mu_cont = arma::accu(postmean.col(j).t() * cor_inv * postmean.col(j));
    var_cont = arma::accu(diag_cor_inv % postvar.col(j));
    obj_val = obj_val - mu_cont / 2.0 - var_cont / 2.0;
  }

  // from penalties

  for (int j = 0; j < nsnps; j++) {
    obj_val = obj_val + pen_bias(bias(j), mean_bias, var_bias);
    obj_val = obj_val + pen_seq_error(seq(j), mean_seq, var_seq);
  }

  return(obj_val);

}
