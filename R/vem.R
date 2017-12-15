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
#' @param seq A vector of initial sequencing errors. Should be
#'     the same length as the number of columns of \code{refmat}
#'     (number of SNPs). Must be between 0 and 1.
#' @param bias A vector of initial bias paramters. Should be the
#'     same length as the number of columns of \code{refmat}
#'     (number of SNPs). Must be greater than 0.
#' @param od A vector of initial overdispersion parameters.
#'     Should be the same length as the number of columns of
#'     \code{refmat} (number of SNPs). Must be between 0 and 1.
#' @param allele_freq A vector of initial allele frequencies.
#'     Should be the same length as the number of columns of
#'     \code{refmat} (number of SNPs). Must be between 0 and 1.
#' @param inbreeding A vector of initial individual-specific
#'     inbreeding coefficients.
#'     Should be the same length as the number of rows of
#'     \code{refmat} (number of individuals).
#' @param cor_mat Initial correlation matrix. Should have the same
#'     number of columns/rows as the number of individuals.
#' @param postmean Initial variational posterior means. Should
#'     have the same dimensions as \code{refmat}.
#' @param postvar Initial posterior variances. Should have the
#'     same dimensions as \code{refmat}.
#' @param update_cor A logical. Should we update the underlying
#'     correlation matrix \code{TRUE} or not \code{FALSE}. Will throw
#'     an error if there are more individuals than SNPs and set to
#'     \code{TRUE}.
#' @param update_inbreeding A logical. Should we update the inbreeding coefficients \code{TRUE}
#'     or not \code{FALSE}?
#' @param verbose Should we print a lot of output \code{TRUE} or not \code{FALSE}?
#' @param control A list of control paramters (\code{itermax}, \code{obj_tol}).
#'
#' @export
#'
#' @author David Gerard
mupdog <- function(refmat,
                   sizemat,
                   ploidy,
                   verbose     = TRUE,
                   mean_bias   = 0,
                   var_bias    = 1,
                   mean_seq    = -4.7,
                   var_seq     = 1,
                   seq         = NULL,
                   bias        = NULL,
                   od          = NULL,
                   allele_freq = NULL,
                   inbreeding  = NULL,
                   cor_mat     = NULL,
                   postmean    = NULL,
                   postvar     = NULL,
                   update_cor  = TRUE,
                   update_inbreeding = TRUE,
                   control = list()) {

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
  if (is.null(seq)) {
    seq <- rep(0.005, length = nsnps)
  }
  if (is.null(bias)) {
    bias <- rep(1, length = nsnps)
  }
  if (is.null(od)) {
    od <- rep(0.001, length = nsnps)
  }
  if (is.null(allele_freq)) {
    allele_freq <- colMeans(refmat / sizemat, na.rm = TRUE)
  }
  if (is.null(inbreeding)) {
    inbreeding <- rep(0, length = nind)
  }
  if (is.null(cor_mat)) {
    cor_mat <- diag(nind)
  }
  if (is.null(postmean)) {
    postmean <- matrix(0, nrow = nind, ncol = nsnps)
  }
  if (is.null(postvar)) {
    postvar <- matrix(1, nrow = nind, ncol = nsnps)
  }
  if (!is.null(control$obj_tol)) {
    obj_tol <- control$obj_tol
  } else {
    obj_tol <- 0.001
  }
  if (!is.null(control$itermax)) {
    itermax <- control$itermax
  } else {
    itermax <- 100
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

  assertthat::assert_that(is.numeric(seq))
  assertthat::are_equal(length(seq), nsnps)
  assertthat::assert_that(all(seq >= 0))
  assertthat::assert_that(all(seq <= 1))

  assertthat::assert_that(is.numeric(bias))
  assertthat::are_equal(length(bias), nsnps)
  assertthat::assert_that(all(bias >= 0))

  assertthat::assert_that(is.numeric(od))
  assertthat::are_equal(length(od), nsnps)
  assertthat::assert_that(all(od >= 0))
  assertthat::assert_that(all(od <= 1))

  assertthat::assert_that(is.numeric(allele_freq))
  assertthat::are_equal(length(allele_freq), nsnps)
  assertthat::assert_that(all(allele_freq >= 0))
  assertthat::assert_that(all(allele_freq <= 1))

  assertthat::assert_that(is.numeric(inbreeding))
  assertthat::are_equal(length(inbreeding), nind)
  assertthat::assert_that(all(inbreeding >= 0))
  assertthat::assert_that(all(inbreeding <= 1))

  assertthat::assert_that(is.numeric(cor_mat))
  assertthat::assert_that(is.matrix(cor_mat))
  assertthat::are_equal(nrow(cor_mat), nind)
  assertthat::are_equal(ncol(cor_mat), nind)
  assertthat::assert_that(all(cor_mat >= 0))
  assertthat::assert_that(all(cor_mat <= 1))
  assertthat::assert_that(all(diag(cor_mat) == 1))
  eigen_decomposition <- eigen(cor_mat, only.values = TRUE)
  assertthat::assert_that(all(eigen_decomposition$values > 0))
  if ((nind > nsnps) & update_cor) {
    stop("if there are more individuals than SNPs, update_core must be set to FALSE.")
  }
  cor_inv <- solve(cor_mat)

  assertthat::assert_that(is.numeric(postmean))
  assertthat::assert_that(is.matrix(postmean))
  assertthat::are_equal(dim(postmean), dim(refmat))

  assertthat::assert_that(is.numeric(postvar))
  assertthat::assert_that(is.matrix(postvar))
  assertthat::are_equal(dim(postvar), dim(refmat))
  assertthat::assert_that(all(postvar > 0))

  assertthat::assert_that(obj_tol > 0)
  assertthat::assert_that(itermax > 0)

  ############################################################

  ## Permanent bounds ------------------------------------------------
  lower_vec <- c(rep(-500, nind), rep(10^-6, nind), 10 ^ -6)
  upper_vec <- c(rep(500, 2 * nind), 1 - 10^-6)

  ## first log-beta array
  lbeta_array <- compute_all_log_bb(refmat = refmat,
                                    sizemat = sizemat,
                                    ploidy = ploidy,
                                    seq = seq,
                                    bias = bias,
                                    od = od)


  obj <- -Inf
  iter <- 1
  err <- Inf
  while (iter < itermax & err > obj_tol) {
    obj_old <- obj
    ## Update posterior means, posterior SD, and allele frequencies -----------------------------------------------
    ## Can parallelize this and also manually calculate gradients for faster computation
    ## Consider using optimx package --- but optim seems faster
    if (verbose) {
      cat(" Updating postmean, postvar, and allele_freq.\n")
      pb <- utils::txtProgressBar(min = 1, max = nsnps, style = 3)
    }
    for (index in 1:nsnps) {
      muSigma2Alpha <- c(postmean[, index], postvar[, index], allele_freq[index])

      oout <- stats::optim(par = muSigma2Alpha, fn = obj_for_mu_wrapper, method = "L-BFGS-B",
                           control = list(fnscale = -1, maxit = 5), lower = lower_vec, upper = upper_vec,
                           rho = inbreeding, log_bb_dense = lbeta_array[, index, ], ploidy = ploidy,
                           cor_inv = cor_inv)

      postmean[, index]  <- oout$par[1:nind]
      postvar[, index]   <- oout$par[(nind + 1):(2 * nind)]
      allele_freq[index] <- oout$par[2 * nind + 1]
      if (verbose) {
        utils::setTxtProgressBar(pb, index)
      }
    }
    if (verbose) {
      close(pb)
    }

    # obj_for_mu(mu = postmean[, index], sigma2 = postvar[, index], alpha = allele_freq[index], rho = inbreeding, log_bb_dense = lbeta_array[, index, ], ploidy = ploidy, cor_inv = cor_inv)
    # obj_for_mu_wrapper(muSigma2Alpha = c(postmean[, index], postvar[, index], allele_freq[index]), rho = inbreeding, log_bb_dense = lbeta_array[, index, ], ploidy = ploidy, cor_inv = cor_inv)

    ## Update inbreeding coefficients ---------------------------------------------------------------------------------------------
    ## Could parallellize this later
    if (update_inbreeding) {
      if (verbose) {
        cat("Done updating postmean, postvar, and allele_freq.\n",
            "Updating inbreeding.\n")
      }
      for (index in 1:nind) {
        oout <- stats::optim(par = inbreeding[index], fn = obj_for_rho, method = "Brent", lower = 0, upper = 1,
                             control = list(fnscale = -1, maxit = 10),
                             mu = postmean[index, ], sigma2 = postvar[index, ], alpha = allele_freq,
                             log_bb_dense = lbeta_array[index, ,], ploidy = ploidy)
        inbreeding[index] <- oout$par
      }
    }


    ## Update correlation matrix ---------------------------------------
    if (update_cor) {
      cor_mat <- update_R(postmean = postmean, postvar = postvar)
      cor_inv <- solve(cor_mat)
    }

    ## Update sequencing error rate, bias, and overdispersion ----------------------------------------------------------------------
    ## can optimize this by actually calculating gradient and parallellizing -------------------------------------------------------
    warray <- compute_all_post_prob(ploidy = ploidy,
                                    mu = postmean,
                                    sigma2 = postvar,
                                    alpha = allele_freq,
                                    rho = inbreeding)

    if (verbose) {
      cat("Done updating inbreeding.\n",
          "Updating seq, bias, and od.\n")
    }
    for (index in 1:nsnps) {
      oout <- stats::optim(par = c(seq[index], bias[index], od[index]), fn = obj_for_eps,
                           method = "L-BFGS-B", lower = rep(10^-6, 3), upper = c(1 - 10^-6, Inf, 1 - 10^-6),
                           control = list(fnscale = -1, maxit = 10),
                           refvec = refmat[, index], sizevec = sizemat[, index], ploidy = ploidy, mean_bias = mean_bias,
                           var_bias = var_bias, mean_seq = mean_seq, var_seq = var_seq, wmat = warray[, index, ])
      seq[index]  <- oout$par[1]
      bias[index] <- oout$par[2]
      od[index]   <- oout$par[3]
    }
    obj_for_eps(parvec = rep(10 ^ -6, 3),
                refvec = refmat[, index], sizevec = sizemat[, index], ploidy = ploidy, mean_bias = mean_bias,
                var_bias = var_bias, mean_seq = mean_seq, var_seq = var_seq, wmat = warray[, index, ])

    lbeta_array <- compute_all_log_bb(refmat = refmat,
                                      sizemat = sizemat,
                                      ploidy = ploidy,
                                      seq = seq,
                                      bias = bias,
                                      od = od)

    if (verbose) {
      cat("Done updating seq, bias, and od.\n\n")
    }

    ## Calculate objective ---------------------------------------------------------------------------------------------------------
    obj <- elbo(warray = warray, lbeta_array = lbeta_array, cor_inv = cor_inv,
                postmean = postmean, postvar = postvar, bias = bias,
                seq = seq, mean_bias = mean_bias, var_bias = var_bias,
                mean_seq = mean_seq, var_seq = var_seq, ploidy = ploidy)

    ## stopping criteria -----------------------------------------------------------------------------------------------------------
    err <- abs(obj_old / obj) - 1

    cat(" iteration:", iter, "\n",
        "objective:", obj, "\n",
        "      err:", err, "\n\n")

    iter <- iter + 1
  }

  map_dosage <- apply(warray, c(1, 2), which.max) - 1

  return_list             <- list()
  return_list$map_dosage  <- map_dosage
  return_list$postprob    <- warray
  return_list$seq         <- seq
  return_list$bias        <- bias
  return_list$od          <- od
  return_list$allele_freq <- allele_freq
  return_list$inbreeding  <- inbreeding
  return_list$cor_mat     <- cor_mat
  return_list$postmean    <- postmean
  return_list$postvar     <- postvar

  return(return_list)
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


