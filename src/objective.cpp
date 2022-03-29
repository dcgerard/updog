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
//' @param sigma2 The variational variance (not standard deviation).
//' @param alpha The allele frequency.
//' @param rho The individual's overdispersion parameter.
//'
//' @return The posterior probability of having \code{dosage} A alleles.
//'
//' @author David
//'
//' @keywords internal
//' @noRd
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
//' @return An array. The rows index the individuals, the columns index the
//'     SNPS, and the third dimension indexes the genotypes. Element (i, j, k)
//'     is the return of \code{\link{post_prob}}.
//'
//' @author David Gerard
//'
//' @keywords internal
//' @noRd
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
//' @keywords internal
//' @noRd
//'
//' @return A three dimensional array. The rows index the individuals, the
//'     columns index the SNPs, and the third dimension indexes the
//'     genotypes. This is the log-likelihood for each individual/snp/genotype
//'     combination.
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


//' Computes \deqn{\Phi^{-1}(F(k|K,\alpha_j,\rho_i))} for all possible (i,j,k).
//'
//' @param alpha A vector whose jth element is the allele frequency of SNP j.
//' @param rho A vector whose ith element is the inbreeding coefficient of individual i.
//' @param ploidy The ploidy of the species.
//'
//' @return A three dimensional array. The rows index the individuals,
//'     the columns index the SNPs, and the third dimension indexes the
//'     genotypes. Computes the "continuous genotype".
//'
//' @author David Gerard
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::Cube<double> compute_all_phifk(NumericVector alpha, NumericVector rho, int ploidy) {
  int nind = rho.length();
  int nsnps = alpha.length();

  arma::Cube<double> phifk(nind, nsnps, ploidy + 2);

  for (int i = 0; i < nind; i++) {
    for (int j = 0; j < nsnps; j++) {
      for (int k = -1; k < ploidy; k++) {
        phifk(i, j, k + 1) = R::qnorm5(pbetabinom_double(k, ploidy, alpha(j), rho(i), false), 0, 1, true, false);
      }
      phifk(i, j, ploidy + 1) = R_PosInf;
    }
  }
  return phifk;
}


//' Penalty on bias parameter.
//'
//' @param h The current value of bias parameter. Must be
//'     greater than 0. A value of 1 means no bias.
//' @param mu_h The prior mean of the log-bias parameter.
//' @param sigma2_h The prior variance (not standard deviation)
//'     of the log-bias parameter. Set to to \code{Inf} to return \code{0}.
//'
//' @return A double. The default penalty on the allelic bias parameter.
//'
//' @author David Gerard
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double pen_bias(double h, double mu_h, double sigma2_h) {

  if (h < TOL) {
    Rcpp::stop("pen_bias: h must be greater than 0.");
  }
  if (sigma2_h < TOL) {
    Rcpp::stop("pen_bias: sigma2_h must be greater tha 0.");
  }

  double pen;
  if (arma::is_finite(sigma2_h)) {
    pen = -std::log(h) - std::pow(std::log(h) - mu_h, 2) / (2.0 * sigma2_h);
  } else {
    pen = 0.0;
  }

  return pen;
}

//' Penalty on sequencing error rate.
//'
//' @param eps The current value of sequencing error rate.
//'     Must be between 0 and 1.
//' @param mu_eps The prior mean of the logit sequencing error rate.
//' @param sigma2_eps The prior variance (not standard deviation)
//'     of the logit sequencing error rate. Set this to \code{Inf} to
//'     return \code{0}.
//'
//' @return A double. The default penalty on the sequencing error rate.
//'
//' @author David Gerard
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double pen_seq_error(double eps, double mu_eps, double sigma2_eps) {

  if ((eps < TOL) | ((1.0 - eps) < TOL)) {
    Rcpp::stop("pen_seq_error: eps must be between 0 and 1.");
  }
  if (sigma2_eps < TOL) {
    Rcpp::stop("pen_seq_error: sigma2_eps must be greater tha 0.");
  }

  double pen;
  if (arma::is_finite(sigma2_eps)) {
    pen = -log(eps * (1.0 - eps)) - std::pow(logit(eps) - mu_eps, 2) / (2.0 * sigma2_eps);
  } else {
    pen = 0.0;
  }
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
//' @return A double. The objective when updating
//'     \code{rho} in \code{\link{mupdog}}.
//'
//'
//' @author David Gerard
//'
//' @keywords internal
//' @noRd
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


//' Objective function when updating alpha
//'
//' @param mu A vector. The ith element is individual i's variational posterior mean at the SNP.
//' @param sigma2 A vector. The ith element is individual i's variational posterior variance at the SNP.
//' @param alpha The SNP's allele frequency.
//' @param rho A vector. The ith element is individuals i's inbreeding coefficient.
//' @param log_bb_dense A matrix of log-densities of the beta binomial. The rows index the individuals and the columns index the allele dosage.
//' @param ploidy The ploidy of the species.
//'
//'
//' @return A double. The objective when updating \code{alpha} in
//'     \code{\link{mupdog}}.
//'
//'
//' @author David Gerard
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double obj_for_alpha(arma::Col<double> mu,
		  arma::Col<double> sigma2,
		  double alpha,
		  NumericVector rho,
		  NumericMatrix log_bb_dense,
		  int ploidy) {

  // check input --------------------------------------------------
  int nind = log_bb_dense.nrow();
  if (log_bb_dense.ncol() != ploidy + 1) {
    Rcpp::stop("obj_for_alpha: log_bb_dense must have ploidy+1 columns.");
  }
  if (mu.n_elem != (unsigned int)nind) {
    Rcpp::Rcout << mu.n_elem << std::endl;
    Rcpp::stop("obj_for_alpha: mu must have length equal to the number of individuals.");
  }
  if (sigma2.n_elem != (unsigned int)nind) {
    Rcpp::Rcout << sigma2.n_elem << std::endl;
    Rcpp::stop("obj_for_alpha: sigma2 must have length equal to the number of individuals.");
  }
  if (rho.length() != nind) {
    Rcpp::stop("obj_for_alpha: rho must have length equal to the number of individuals.");
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
//' @param mean_od The prior mean of the logit of the overdispersion parameter
//' @param var_od The prior variance of the logit of the overdispersion parameter.
//' @param wmat The matrix of (variational) posterior probabilities for each dosage.
//'     The rows index the individuals and the columns index the dosage levels.
//' @param update_seq A logical. This is not used in \code{obj_for_eps},
//'     but sets the first element to \code{0.0} in \code{\link{grad_for_eps}}.
//' @param update_bias A logical. This is not used in \code{obj_for_eps},
//'     but sets the second element to \code{0.0} in \code{\link{grad_for_eps}}.
//' @param update_od A logical. This is not used in \code{obj_for_eps},
//'     but sets the third element to \code{0.0} in \code{\link{grad_for_eps}}.
//'
//' @return A double. The objective when updating \code{eps} in
//'     \code{\link{mupdog}}.
//'
//' @author David Gerard
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double obj_for_eps(NumericVector parvec,
		   NumericVector refvec,
		   NumericVector sizevec,
		   int ploidy,
		   double mean_bias,
		   double var_bias,
		   double mean_seq,
		   double var_seq,
		   double mean_od,
		   double var_od,
		   NumericMatrix wmat,
		   bool update_bias = true,
		   bool update_seq = true,
		   bool update_od = true) {
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
      if (!R_IsNA(refvec(i)) && !R_IsNA(sizevec(i))) {
	xi_current = xi_double((double)k / (double)ploidy, eps, h);
	obj_val = obj_val + wmat(i, k) * dbetabinom_double(refvec(i), sizevec(i), xi_current, tau, true);
      }
    }
  }

  obj_val = obj_val +
    pen_bias(h, mean_bias, var_bias) +
    pen_seq_error(eps, mean_seq, var_seq) +
    pen_seq_error(tau, mean_od, var_od); // can use pen_seq_error() for overdispersion as well!

  if (obj_val == R_NegInf) {
    Rcpp::Rcout << obj_val << std::endl;
    Rcpp::Rcout << parvec << std::endl;
  }

  return obj_val;
}


//' Objective function when updating mu and sigma2.
//'
//' @param mu A vector, the ith element is the variational posterior mean of individual i for the SNP.
//' @param sigma2 A vector, the ith element is the variational posterior variance of individual i for the SNP.
//' @param phifk_mat A matrix that contains the standard normal quantile of the beta-binomial cdf at dosage k for individual i.
//'     The rows index the individuals and the columns index the dosages.
//' @param cor_inv The inverse of the underlying correlation matrix.
//' @param log_bb_dense A matrix of log-densities of the beta binomial. The rows index the individuals and the columns index the allele dosage.
//'     Allele dosage goes from -1 to ploidy, so there are ploidy + 2 elements.
//'
//' @return A double. The objective when updating \code{mu} and \code{sigma2}
//'     in \code{\link{mupdog}}.
//'
//' @author David Gerard
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double obj_for_mu_sigma2(arma::Col<double> mu, arma::Col<double> sigma2, NumericMatrix phifk_mat,
                         arma::Mat<double> cor_inv, NumericMatrix log_bb_dense) {
  // Check input --------------------------------------------------------------------------------
  int nind = log_bb_dense.nrow();
  int ploidy = log_bb_dense.ncol() - 1;
  if (mu.n_elem != (unsigned int)nind) {
    Rcpp::stop("obj_for_mu_sigma2: mu needs to have the same length as the number of columns in log_bb_dense.");
  }
  if (sigma2.n_elem != (unsigned int)nind) {
    Rcpp::stop("obj_for_mu_sigma2: sigma2 needs to have the same length as the number of columns in log_bb_dense.");
  }
  if (phifk_mat.nrow() != nind) {
    Rcpp::stop("obj_for_mu_sigma2: phifk_mat and log_bb_dense needs to have the same dimensions.");
  }
  if (phifk_mat.ncol() != ploidy + 2) {
    Rcpp::stop("obj_for_mu_sigma2: phifk_mat needs to have ploidy+2 columns.");
  }
  if (cor_inv.n_cols != (unsigned int)nind) {
    Rcpp::stop("obj_for_mu_sigma2: cor_inv needs to have the same number of rows as log_bb_dense.");
  }
  if (cor_inv.n_rows != (unsigned int)nind) {
    Rcpp::stop("obj_for_mu_sigma2: cor_inv needs to have the same number of columns as rows.");
  }

  double obj_val = 0.0;
  double current_weight = 0.0;
  for (int i = 0; i < nind; i++) {
    for (int k = 0; k <= ploidy; k++) {
      if (!R_IsNA(log_bb_dense(i, k))) {
        current_weight = R::pnorm5((phifk_mat(i, k + 1) - mu(i)) / std::sqrt(sigma2(i)), 0.0, 1.0, true, false) -
          R::pnorm5((phifk_mat(i, k) - mu(i)) / std::sqrt(sigma2(i)), 0.0, 1.0, true, false);
        obj_val = obj_val + current_weight * log_bb_dense(i, k);
      }
    }
  }

  // Multivariate normal prior integral -------------------------------------
  double muRmu = arma::accu(mu.t() * cor_inv * mu);
  double trRsigma2 = arma::accu(arma::diagvec(cor_inv).t() * sigma2);
  obj_val = obj_val - muRmu / 2.0 - trRsigma2 / 2.0 + arma::accu(arma::log(sigma2)) / 2.0;

  return obj_val;
}

//' Wrapper for \code{\link{obj_for_mu_sigma2}} so that I can use it in \code{optim}.
//'
//' @inheritParams obj_for_mu_sigma2
//' @param muSigma2 A vector. The first half are mu and the second half are sigma2.
//' @inherit obj_for_mu_sigma2 return
//'
//'
//' @author David Gerard
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double obj_for_mu_sigma2_wrapper(arma::Col<double> muSigma2, NumericMatrix phifk_mat,
                                 arma::Mat<double> cor_inv, NumericMatrix log_bb_dense) {
  int nind = muSigma2.n_elem / 2;
  double obj = obj_for_mu_sigma2(muSigma2.head_rows(nind), muSigma2.tail_rows(nind), phifk_mat,
                                 cor_inv, log_bb_dense);
  return obj;
}


//' The evidence lower bound
//'
//' @param warray An three-way array. The (i,j,k)th entry is the variational posterior probability
//'     that individual i at SNP j has dosage k - 1. See \code{\link{compute_all_post_prob}}.
//' @param lbeta_array A three-way array. The (i,j,k)th entry is the log-density of the betabinomial
//'     for individual i at SNP j and dosage k - 1. See \code{\link{compute_all_log_bb}}.
//' @param cor_inv The inverse of the correlation matrix.
//' @param postmean A matrix. The (i,j)th entry is the variational posterior mean for individual i
//'     at SNP j.
//' @param postvar A matrix. The (i,j)th entry is the variational posterior variance for individual i at SNP j.
//' @param bias A vector. The jth entry is the allele bias for SNP j.
//' @param seq A vector. The jth entry is the sequencing error rate at SNP j.
//' @param mean_bias The prior mean on the log-bias.
//' @param var_bias The prior variance on the log-bias.
//' @param mean_seq The prior mean on the logit of the sequencing error rate.
//' @param var_seq The prior variance on the logit of the sequencing error rate.
//' @param ploidy The ploidy of the species.
//'
//' @return A double. The evidence lower-bound that \code{\link{mupdog}}
//'     maximizes.
//'
//'
//' @author David Gerard
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double elbo(arma::Cube<double> warray,
            arma::Cube<double> lbeta_array,
            arma::Mat<double> cor_inv,
            arma::Mat<double> postmean,
            arma::Mat<double> postvar,
            NumericVector bias,
            NumericVector seq,
            double mean_bias,
            double var_bias,
            double mean_seq,
            double var_seq,
            int ploidy) {
  // Check input -------------------------------------------------------
  int nsnps = warray.n_cols;
  int nind  = warray.n_rows;
  if (lbeta_array.n_cols != (unsigned int)nsnps) {
    Rcpp::stop("elbo: lbeta_array must have the same dimension as warray");
  }
  if (lbeta_array.n_rows != (unsigned int)nind) {
    Rcpp::stop("elbo: lbeta_array must have the same dimension as warray");
  }
  if (lbeta_array.n_slices != ((unsigned int)ploidy + 1)) {
    Rcpp::stop("elbo: lbeta_array must have third-dimension ploidy+1.");
  }
  if (warray.n_slices != ((unsigned int)ploidy + 1)) {
    Rcpp::stop("elbo: warray must have third-dimension ploidy+1.");
  }
  if (cor_inv.n_rows != (unsigned int)nind) {
    Rcpp::stop("elbo: cor_inv must have the number of rows equal to the number of individuals.");
  }
  if (cor_inv.n_cols != (unsigned int)nind) {
    Rcpp::stop("elbo: cor_inv must have the number of columns equal to the number of individuals.");
  }
  if (postmean.n_rows != (unsigned int)nind) {
    Rcpp::stop("elbo: postmean must have the number of rows equal to the number of individuals.");
  }
  if (postmean.n_cols != (unsigned int)nsnps) {
    Rcpp::stop("elbo: postmean must have the number of columns equal to the number of SNPs.");
  }
  if (postvar.n_rows != (unsigned int)nind) {
    Rcpp::stop("elbo: postvar must have the number of rows equal to the number of individuals.");
  }
  if (postvar.n_cols != (unsigned int)nsnps) {
    Rcpp::stop("elbo: postvar must have the number of columns equal to the number of SNPs.");
  }
  if (bias.length() != nsnps) {
    Rcpp::stop("elbo: bias must have number of elements equal to the number of SNPs.");
  }
  if (seq.length() != nsnps) {
    Rcpp::stop("elbo: seq must have number of elements equal to the number of SNPs.");
  }
  if (var_bias < 0) {
    Rcpp::stop("elbo: var_bais must be greater than 0.");
  }
  if (var_seq < 0) {
    Rcpp::stop("elbo: var_seq must be greater than 0.");
  }

  // Objective Function ---------------------------------------------------------------

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



//' Objective function for updating the beta-binomial genotype distribution when
//' \code{model = "bb"} in \code{\link{flex_update_pivec}}.
//'
//' @param parvec A vector of length 2. The first term is the current mean of the
//'     underlying beta. The second term is the current overdispersion parameter.
//' @param ploidy The ploidy of the species.
//' @param weight_vec A vector of length \code{ploidy + 1} that contains the weights
//'     for each component beta-binomial.
//'
//' @return A double. The objective when updating the beta-binomial genotype
//'     distribution in \code{\link{mupdog}}.
//'
//' @author David Gerard
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double obj_for_weighted_lbb(NumericVector parvec,
                            int ploidy,
                            NumericVector weight_vec) {
  if (parvec.length() != 2) {
    Rcpp::stop("obj_for_weighted_lbb: parvec not of length 2.");
  }
  if (weight_vec.length() != (ploidy + 1)) {
    Rcpp::stop("obj_for_weighted_lbb: weight_vec not of length ploidy + 1.");
  }

  double mu  = parvec(0);
  double tau = parvec(1);
  double obj = 0.0;
  for (int i = 0; i <= ploidy; i++) {
    obj = obj + weight_vec(i) * dbetabinom_double(i, ploidy, mu, tau, true);
  }
  return obj;
}


//' Objective function for updating discrete normal genotype distribution
//' when \code{model = "normal"} in \code{\link{flex_update_pivec}}.
//'
//' @param parvec A vector of length 2. The first term is the current mean of the
//'     underlying normal. The second term is the current standard deviation
//'     (not variance) of the normal.
//' @param ploidy The ploidy of the species.
//' @param weight_vec A vector of length \code{ploidy + 1} that contains the weights
//'     for each component beta-binomial.
//'
//' @return A double. The objective when updating the normal
//'     genotype distribution in \code{\link{mupdog}}.
//'
//' @author David Gerard
//'
//' @keywords internal
//' @noRd
//'
// [[Rcpp::export]]
double obj_for_weighted_lnorm(NumericVector parvec,
                              int ploidy,
                              NumericVector weight_vec) {
  if (parvec.length() != 2) {
    Rcpp::stop("obj_for_weighted_lbb: parvec not of length 2.");
  }
  if (weight_vec.length() != (ploidy + 1)) {
    Rcpp::stop("obj_for_weighted_lbb: weight_vec not of length ploidy + 1.");
  }

  double mu    = parvec(0);
  double sigma = parvec(1);
  NumericVector lpvec(ploidy + 1);
  double obj = 0.0;
  for (int i = 0; i <= ploidy; i++) {
    lpvec(i) = R::dnorm((double)i, mu, sigma, true);
    obj = obj + weight_vec(i) * lpvec(i);
  }
  double lsum = log_sum_exp(lpvec);
  obj = obj - Rcpp::sum(weight_vec) * lsum;
  return(obj);
}







