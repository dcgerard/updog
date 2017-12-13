### Code for VEM

#' Multivariate updog
#'
#' @param refmat A matrix of reference counts.
#'     The rows index the individuals and the columns index the SNPs.
#' @param sizemat A matrix of total counts.
#'     The rows index the individuals and the columns index the SNPs.
#' @param ploidy The ploidy of the species. To estimate the ploidy,
#'     re-run updog at various ploidy levels and choose the one with
#'     the largest ELBO. This assumes that the ploidy is the same for
#'     all individuals in the sample.
#' @param mean_bias The prior mean of the log-bias. Defaults to 0
#'     (no bias).
#' @param var_bias The prior variance on the log-bias. Defaults to 1.
#'     This roughly corresponds to likely bias values between 0.14 and
#'     7.4. This is a far wider interval than what we observe in
#'     practice, thus making this prior rather uniformative. We usually
#'     observe bias values somewhere between 0.5 and 2.
#' @param mean_seq The prior mean of the logit-sequencing-error-rate.
#'     Defaults to -4.7. This corresponds to a sequencing error rate of
#'     0.009.
#' @param var_seq The prior variance of the logit-sequencing-error-rate.
#'     Defaults to 1. This corresponds to likely values of 0.001 and
#'     0.06. This upper bound is larger than what we would expect given
#'     the current state of next-gen-sequencings.
#' @param seq_init A vector of initial sequencing errors. Should be
#'     the same length as the number of columns of \code{refmat}
#'     (number of SNPs). Must be between 0 and 1.
#' @param bias_init A vector of initial bias paramters. Should be the
#'     same length as the number of columns of \code{refmat}
#'     (number of SNPs). Must be greater than 0.
#' @param od_init A vector of initial overdispersion parameters.
#'     Should be the same length as the number of columns of
#'     \code{refmat} (number of SNPs). Must be between 0 and 1.
#' @param allele_freq_init A vector of initial allele frequencies.
#'     Should be the same length as the number of columns of
#'     \code{refmat} (number of SNPs). Must be between 0 and 1.
#' @param inbreeding_init A vector of initial individual-specific
#'     inbreeding coefficients.
#'     Should be the same length as the number of rows of
#'     \code{refmat} (number of individuals).
#' @param cor_init Initial correlation matrix. Should have the same
#'     number of columns/rows as the number of individuals.
#' @param postmean_init Initial variational posterior means. Should
#'     have the same dimensions as \code{refmat}.
#' @param postvar_init Initial posterior variances. Should have the
#'     same dimensions as \code{refmat}.
#' @param update_cor A logical. Should we update the underlying
#'     correlation matrix \code{TRUE} or not \code{FALSE}. Will throw
#'     an error if there are more individuals than SNPs and set to
#'     \code{TRUE}.
#'
#'
#' @author David Gerard
mupdog <- function(refmat,
                   sizemat,
                   ploidy,
                   mean_bias        = 0,
                   var_bias         = 1,
                   mean_seq         = -4.7,
                   var_seq          = 1,
                   seq_init         = NULL,
                   bias_init        = NULL,
                   od_init          = NULL,
                   allele_freq_init = NULL,
                   inbreeding_init  = NULL,
                   cor_init         = NULL,
                   postmean_init    = NULL,
                   postvar_init     = NULL,
                   update_cor       = TRUE) {

  ##########################################################
  ## Check refmat and sizemat ------------------------------
  ##########################################################
  assertthat::assert_that(is.matrix(refmat))
  assertthat::assert_that(is.numeric(refmat))
  assertthat::assert_that(is.matrix(sizemat))
  assertthat::assert_that(is.numeric(sizemat))
  assertthat::are_equal(dim(refmat), dim(sizemat))
  assertthat::assert_that(all(sizemat >= refmat, na.rm = TRUE))
  assertthat::assert_that(all(sizemat >= 0, na.rm = TRUE))

  nind <- nrow(refmat)
  nsnps <- ncol(refmat)

  ###########################################################
  ## Setup Default parameters -------------------------------
  ###########################################################
  if (is.null(seq_init)) {
    seq_init <- rep(0.005, length = nsnps)
  }
  if (is.null(bias_init)) {
    bias_init <- rep(0, length = nsnps)
  }
  if (is.null(od_init)) {
    od_init <- rep(0.001, length = nsnps)
  }
  if (is.null(allele_freq_init)) {
    allele_freq_init <- colMeans(refmat / sizemat, na.rm = TRUE)
  }
  if (is.null(inbreeding_init)) {
    inbreeding_init <- rep(0, length = nind)
  }
  if (is.null(cor_init)) {
    cor_init <- diag(nind)
  }
  if (is.null(postmean_init)) {
    postmean_init <- matrix(0, nrow = nind, ncol = nsnps)
  }
  if (is.null(postvar_init)) {
    postvar_init <- matrix(1, nrow = nind, ncol = nsnps)
  }

  ######################################################
  ## Check rest of input -------------------------------
  ######################################################

  assertthat::assert_that(is.numeric(mean_bias))
  assertthat::are_equal(length(mean_bias), 1)

  assertthat::assert_that(is.numeric(var_bias))
  assertthat::are_equal(length(var_bias), 1)
  assertthat::assert_that(var_bias > 0)

  assertthat::assert_that(is.numeric(mean_seq))
  assertthat::are_equal(length(mean_seq), 1)

  assertthat::assert_that(is.numeric(var_seq))
  assertthat::are_equal(length(var_seq), 1)
  assertthat::assert_that(var_seq > 0)

  assertthat::assert_that(is.numeric(seq_init))
  assertthat::are_equal(length(seq_init), nsnps)
  assertthat::assert_that(all(seq_init >= 0))
  assertthat::assert_that(all(seq_init <= 1))

  assertthat::assert_that(is.numeric(bias_init))
  assertthat::are_equal(length(bias_init), nsnps)
  assertthat::assert_that(all(bias_init >= 0))

  assertthat::assert_that(is.numeric(od_init))
  assertthat::are_equal(length(od_init), nsnps)
  assertthat::assert_that(all(od_init >= 0))
  assertthat::assert_that(all(od_init <= 1))

  assertthat::assert_that(is.numeric(allele_freq_init))
  assertthat::are_equal(length(allele_freq_init), nsnps)
  assertthat::assert_that(all(allele_freq_init >= 0))
  assertthat::assert_that(all(allele_freq_init <= 1))

  assertthat::assert_that(is.numeric(inbreeding_init))
  assertthat::are_equal(length(inbreeding_init), nind)
  assertthat::assert_that(all(inbreeding_init >= 0))
  assertthat::assert_that(all(inbreeding_init <= 1))

  assertthat::assert_that(is.numeric(cor_init))
  assertthat::assert_that(is.matrix(cor_init))
  assertthat::are_equal(nrow(cor_init), nind)
  assertthat::are_equal(ncol(cor_init), nind)
  assertthat::assert_that(all(cor_init >= 0))
  assertthat::assert_that(all(cor_init <= 1))
  assertthat::assert_that(all(diag(cor_init) == 1))
  eigen_decomposition <- eigen(cor_init, only.values = TRUE)
  assertthat::assert_that(all(eigen_decomposition$values > 0))
  if ((nind > nsnps) & update_cor) {
    stop("if there are more individuals than SNPs, update_core must be set to FALSE.")
  }

  assertthat::assert_that(is.numeric(postmean_init))
  assertthat::assert_that(is.matrix(postmean_init))
  assertthat::are_equal(dim(postmean_init), dim(refmat))

  assertthat::assert_that(is.numeric(postvar_init))
  assertthat::assert_that(is.matrix(postvar_init))
  assertthat::are_equal(dim(postvar_init), dim(refmat))
  assertthat::assert_that(all(postvar_init > 0))

  ############################################################

  postmean    <- postmean_init
  postvar     <- postvar_init
  seq_error   <- seq_init
  bias        <- bias_init
  od          <- od_init
  allele_freq <- allele_freq_init
  inbreeding  <- inbreeding_init
  cor_mat     <- cor_init
  postmean    <- postmean_init
  postvar     <- postvar_init


  if (update_cor) {
    cor_mat <- update_R(postmean = postmean, postvar = postvar)
  }



}

#' Update the underlying correlation matrix.
#'
#' @param postmean The matrix of posterior means. The rows index the
#'     individuals and the columns index the SNPs.
#' @param postvar The matrix of posterior variances. The rows index
#'     the individuals and the columns index the SNPs.
#'
#' @author David Gerard
update_R <- function(postmean, postvar) {
  stats::cov2cor(tcrossprod(postmean) + diag(rowSums(postvar)))
}


