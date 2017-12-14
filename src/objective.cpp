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
    Rcpp::stop("sizemat and refmat must have the same dimensions.");
  }
  if (sizemat.ncol() != nsnps) {
    Rcpp::stop("sizemat and refmat must have the same dimensions.");
  }
  if (seq.length() != nsnps) {
    Rcpp::stop("seq must have the same length as the number of SNPs.");
  }
  if (bias.length() != nsnps) {
    Rcpp::stop("bias must have the same length as the number of SNPs.");
  }
  if (od.length() != nsnps) {
    Rcpp::stop("od must have the same length as the number of SNPs.");
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
    Rcpp::stop("log_bb_dense must have ploidy+1 columns.");
  }
  if (mu.length() != nsnps) {
    Rcpp::stop("mu must have length equal to the number of SNPs");
  }
  if (sigma2.length() != nsnps) {
    Rcpp::stop("sigma2 must have length equal to the number of SNPs");
  }
  if (alpha.length() != nsnps) {
    Rcpp::stop("alpha must have length equal to the number of SNPs");
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
