### Code for VEM

#' Multi-SNP updog.
#'
#' A method to genotype autopolyploids using GBS or RAD-seq like data by accounting
#' for correlations in the genotype distribution between the individuals.
#'
#' Blischak et al (2017) developed a genotyping approach for autopolyploids
#' that assumes a Balding-Nichols generative model (Balding and Nichols, 1997)
#' on the genotypes. Using a different generative model, Gerard et al (2018)
#' accounted for common issues in sequencing data ignored by previous researchers.
#' Mupdog unites and extends these two approaches:
#' \itemize{
#' \item{Unite: We account for locus-specific allele-bias, locus-specific sequencing error,
#'       and locus-specific overdispersion while marginally assuming a Balding-Nichols generative model on the genotypes.}
#' \item{Extend: We account for underlying correlations between the individuals using a Gaussian copula model.}
#' }
#' Mupdog uses a variational Bayes approach to estimate all parameters of interest and
#' the posterior probabilities of the genotypes for each individual at each locus.
#'
#' @param refmat A matrix of reference counts.
#'     The rows index the individuals and the columns index the SNPs.
#' @param sizemat A matrix of total counts.
#'     The rows index the individuals and the columns index the SNPs.
#'     Should have the same dimensions as \code{refmat}.
#' @param ploidy The ploidy of the species. To estimate the ploidy,
#'     re-run \code{mupdog} at various ploidy levels and choose the one with
#'     the largest ELBO. This assumes that the ploidy is the same for
#'     all individuals in the sample.
#' @param mean_bias The prior mean of the log-bias. Defaults to 0
#'     (no bias).
#' @param var_bias The prior variance on the log-bias. Defaults to 1.
#'     This roughly corresponds to likely bias values between 0.14 and
#'     7.4. This is a far wider interval than what we observe in
#'     practice, thus making this prior rather uninformative. We usually
#'     observe bias values somewhere between 0.5 and 2.
#' @param mean_seq The prior mean of the logit-sequencing-error-rate.
#'     Defaults to -4.7. This corresponds to a sequencing error rate of
#'     0.009.
#' @param var_seq The prior variance of the logit-sequencing-error-rate.
#'     Defaults to 1. This corresponds to likely values of 0.001 and
#'     0.06. This upper bound is larger than what we would expect given
#'     the current state of next-gen-sequencing technology.
#' @param seq A vector of initial sequencing errors. Should be
#'     the same length as the number of columns of \code{refmat}
#'     (number of SNPs). Must be between 0 and 1.
#' @param bias A vector of initial bias parameters. Should be the
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
#'     \code{refmat} (number of individuals). Must be between 0 and 1.
#' @param cor_mat Initial correlation matrix. Should have the same
#'     number of columns/rows as the number of individuals.
#' @param postmean Initial variational posterior means. Should
#'     have the same dimensions as \code{refmat}.
#' @param postvar Initial posterior variances. Should have the
#'     same dimensions as \code{refmat}.
#' @param update_cor A logical. Should we update the underlying
#'     correlation matrix \code{TRUE} or not \code{FALSE}. It might be
#'     unwise to set this to \code{TRUE} if you have more individuals
#'     than SNPs.
#' @param update_inbreeding A logical. Should we update the
#'     inbreeding coefficients \code{TRUE}
#'     or not \code{FALSE}?
#' @param update_allele_freq A logical. Should we update the
#'     allele frequencies \code{TRUE} or not \code{FALSE}?
#' @param verbose Should we print a lot of output \code{TRUE}
#'     or not \code{FALSE}?
#' @param control A list of control parameters (\code{itermax},
#'     \code{obj_tol}).
#' @param update_method What generic optimizer should we use to update \code{allele_freq}
#'     and \code{inbreeding}? Options are either \code{"Brent"} or \code{"L-BFGS-B"}.
#'     See \code{\link[stats]{optim}} for details on these optimizers.
#' @param num_core The number of cores to use if you want to
#'     run the optimization steps in parallel. If \code{num_core = 1},
#'     then the optimization step will not be run in parallel.
#'
#' @references
#' \itemize{
#'   \item{David J Balding and Richard A Nichols. \href{http://dx.doi.org/10.1038/sj.hdy.6881750}{Significant genetic correlations among caucasians at forensic DNA loci}. \emph{Heredity}, 78(6):583–589, 1997. doi: 10.1038/sj.hdy.6881750.}
#'   \item{Paul D Blischak, Laura S Kubatko, and Andrea D Wolfe. \href{http://dx.doi.org/10.1093/bioinformatics/btx587}{SNP genotyping and parameter estimation in polyploids using low-coverage sequencing data}. \emph{Bioinformatics}, page btx587, 2017. doi: 10.1093/bioinformatics/btx587.}
#'   \item{Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018). Genotyping Polyploids from Messy Sequencing Data. \emph{Genetics}, 210(3), 789-807. doi: \href{https://doi.org/10.1534/genetics.118.301468}{10.1534/genetics.118.301468}.}
#' }
#'
#' @export
#'
#' @author David Gerard
#'
#' @return A list with some or all of the following elements:
#' \describe{
#'     \item{\code{map_dosage}}{A matrix of numerics containing
#'         the variational maximum-a-posterior (MAP) genotypes
#'         (allele dosages) for each individual at each SNP. Element
#'         (i, j) is the MAP genotype for individual i at SNP j.}
#'     \item{\code{maxpostprob}}{A matrix of numerics containing
#'         the variational maximum posterior probabilities for each
#'         individual at each SNP. The (i, j)th element is the
#'         variational posterior
#'         probability that individual i is genotyped correctly at
#'         SNP j.}
#'     \item{\code{postprob}}{A three-way array of numerics. Element (i, j, k) is
#'         the variational posterior probability that individual i has genotype
#'         k-1 at SNP j.}
#'     \item{\code{seq}}{A vector of numerics. Element j is the estimated
#'         sequencing error rate for SNP j.}
#'     \item{\code{bias}}{A vector of numerics. Element j is the estimated
#'         allelic bias for SNP j.}
#'     \item{\code{od}}{A vector of numerics. Element j is the estimated
#'         overdispersion parameter for SNP j.}
#'     \item{\code{allele_freq}}{A vector of numerics. Element j is the
#'         estimated allele-frequency for SNP j.}
#'     \item{\code{inbreeding}}{A vector of numerics. Element i is the
#'         estimated inbreeding coefficient for individual i.}
#'     \item{\code{cor_mat}}{A symmetric matrix of numerics. Element (i, j)
#'         is the estimated _latent_ correlation between individual
#'         i and individual j.}
#'     \item{\code{postmean}}{A matrix of numerics. Element (i, j) is the
#'         variational posterior mean for individual i at SNP j.}
#'     \item{\code{postvar}}{A matrix of numerics. Element (i, j) is the
#'         variational posterior variance for individual i at SNP j.}
#'     \item{\code{input$refmat}}{A matrix of numerics.
#'         The inputted \code{refmat}.}
#'     \item{\code{input$sizemat}}{A matrix of numerics. The inputted
#'         \code{sizemat}.}
#'     \item{\code{input$ploidy}}{The inputted \code{ploidy}.}
#'     \item{\code{obj}}{The maximized variational objective.}
#' }
#'
#'
#' @examples
#' \donttest{
#' data(uitdewilligen)
#' mout <- mupdog(refmat = uitdewilligen$refmat,
#'                sizemat = uitdewilligen$sizemat,
#'                ploidy = uitdewilligen$ploidy,
#'                verbose = FALSE,
#'                control = list(obj_tol = 10^-4))
#'
#' ## Summaries of output
#' plot(mout, 4)
#' hist(mout$bias)
#' hist(mout$seq)
#' hist(mout$od)
#' hist(mout$inbreeding)
#' hist(mout$allele_freq)
#'
#' ## mupdog can correctly estimate ploidy to be 4
#' mout2 <- mupdog(refmat = uitdewilligen$refmat,
#'                 sizemat = uitdewilligen$sizemat,
#'                 ploidy = 2,
#'                 verbose = FALSE,
#'                 control = list(obj_tol = 10^-4))
#'
#' mout6 <- mupdog(refmat = uitdewilligen$refmat,
#'                 sizemat = uitdewilligen$sizemat,
#'                 ploidy = 6,
#'                 verbose = FALSE,
#'                 control = list(obj_tol = 10^-4))
#'
#' mout8 <- mupdog(refmat = uitdewilligen$refmat,
#'                 sizemat = uitdewilligen$sizemat,
#'                 ploidy = 8,
#'                 verbose = FALSE,
#'                 control = list(obj_tol = 10^-4))
#'
#' y <- c(mout2$obj, mout$obj, mout6$obj, mout8$obj)
#' x <- seq(2, 8, by = 2)
#' plot(x, y, type = "l", xlab = "ploidy", ylab = "objective")
#' }
#'
#'
mupdog <- function(refmat,
                   sizemat,
                   ploidy,
                   verbose            = TRUE,
                   mean_bias          = 0,
                   var_bias           = 1,
                   mean_seq           = -4.7,
                   var_seq            = 1,
                   seq                = NULL,
                   bias               = NULL,
                   od                 = NULL,
                   allele_freq        = NULL,
                   inbreeding         = NULL,
                   cor_mat            = NULL,
                   postmean           = NULL,
                   postvar            = NULL,
                   update_cor         = TRUE,
                   update_inbreeding  = TRUE,
                   update_allele_freq = TRUE,
                   num_core           = 1,
                   update_method      = c("Brent", "L-BFGS-B"),
                   control            = list()) {

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

  update_method <- match.arg(update_method)

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
  eigen_decomposition <- eigen(cor_mat, only.values = TRUE)
  assertthat::assert_that(all(eigen_decomposition$values > 0))
  cor_inv <- solve(cor_mat)

  assertthat::assert_that(is.numeric(postmean))
  assertthat::assert_that(is.matrix(postmean))
  assertthat::are_equal(dim(postmean), dim(refmat))

  assertthat::assert_that(is.numeric(postvar))
  assertthat::assert_that(is.matrix(postvar))
  assertthat::are_equal(dim(postvar), dim(refmat))
  assertthat::assert_that(all(postvar > 0))

  assertthat::assert_that(obj_tol > 0)
  assertthat::assert_that(itermax > 1)
  assertthat::assert_that(num_core >= 1)
  assertthat::assert_that(is.logical(update_allele_freq))
  assertthat::assert_that(is.logical(update_cor))
  assertthat::assert_that(is.logical(update_inbreeding))

  ############################################################

  ## Permanent bounds ------------------------------------------------
  lower_vec <- c(rep(-500, nind), rep(10^-6, nind))
  upper_vec <- c(rep(500, 2 * nind))

  ## first log-beta array
  lbeta_array <- compute_all_log_bb(refmat = refmat,
                                    sizemat = sizemat,
                                    ploidy = ploidy,
                                    seq = seq,
                                    bias = bias,
                                    od = od)

  phifk_array <- compute_all_phifk(alpha = allele_freq,
                                   rho = inbreeding,
                                   ploidy = ploidy)


  obj <- -Inf
  iter <- 1
  err <- Inf
  if (num_core == 1) { # only need to do this once if no cores.
    foreach::registerDoSEQ()
  } else if (num_core > 1) {
    cl = parallel::makeCluster(num_core)
    doParallel::registerDoParallel(cl = cl)
    if (foreach::getDoParWorkers() == 1) {
      stop("num_core > 1 but only one core registered using foreach::registerDoSEQ.")
    }
  } else {
    stop("mupdog: how did you get here?")
  }
  while (iter < itermax & err > obj_tol) {
    obj_old <- obj


    ## Update variational means and variances ---------------------------------
    if (verbose) {
      cat("Updating posterior means and variances.\n\n")
    }

    index <- NA ## stupid workaround to get rid of CRAN note.
    fout <- foreach::foreach(index = seq_len(nsnps), .combine = cbind,
                             .export = c("obj_for_mu_sigma2_wrapper",
                                         "grad_for_mu_sigma2_wrapper")) %dopar% {
      oout <- stats::optim(par = c(postmean[, index], postvar[, index]),
                           fn = obj_for_mu_sigma2_wrapper,
                           gr = grad_for_mu_sigma2_wrapper,
                           method = "L-BFGS-B",
                           lower = lower_vec,
                           upper = upper_vec,
                           control = list(fnscale = -1, maxit = 10),
                           phifk_mat = phifk_array[, index, ],
                           cor_inv = cor_inv,
                           log_bb_dense = lbeta_array[, index, ])
      oout$par
    }
    postmean <- fout[seq_len(nind), ]
    postvar  <- fout[(nind + 1):(2 * nind), ]


    ## Update allele frequencies -----------------------------------------------

    if (update_allele_freq) {
      if (verbose) {
        cat("Updating allele_freq.\n\n")
      }
      allele_freq <- foreach::foreach(index = seq_len(nsnps), .combine = c,
                       .export = "obj_for_alpha") %dopar% {
        oout <- stats::optim(par = allele_freq[index],
                             fn = obj_for_alpha,
                             method = update_method,
                             control = list(fnscale = -1, maxit = 10),
                             lower = 0, upper = 1,
                             mu = postmean[, index],
                             sigma2 = postvar[, index],
                             rho = inbreeding,
                             log_bb_dense = lbeta_array[, index, ],
                             ploidy = ploidy)
        oout$par
        }
    }



    ## Update inbreeding coefficients ---------------------------------------------------------------------------------------------
    if (update_inbreeding) {
      if (verbose) {
        cat("Updating inbreeding.\n\n")
      }
      inbreeding <- foreach::foreach(index = seq_len(nind), .combine = c,
                                     .export = "obj_for_rho") %dopar% {
        oout <- stats::optim(par = inbreeding[index],
                             fn = obj_for_rho,
                             method = update_method,
                             lower = 0, upper = 1,
                             control = list(fnscale = -1, maxit = 10),
                             mu = postmean[index, ],
                             sigma2 = postvar[index, ],
                             alpha = allele_freq,
                             log_bb_dense = lbeta_array[index, ,],
                             ploidy = ploidy)
        oout$par
      }
    }

    ## Recalculate phifk_array after updating allele_freq and inbreeding ---------------------------------------------------------
    phifk_array <- compute_all_phifk(alpha = allele_freq,
                                     rho = inbreeding,
                                     ploidy = ploidy)


    ## Update correlation matrix ---------------------------------------
    if (update_cor) {
      cor_mat <- update_R(postmean = postmean, postvar = postvar)
      cor_inv <- solve(cor_mat)
    }

    ## Update sequencing error rate, bias, and overdispersion ----------------------------------------------------------------------
    warray <- compute_all_post_prob(ploidy = ploidy,
                                    mu = postmean,
                                    sigma2 = postvar,
                                    alpha = allele_freq,
                                    rho = inbreeding)

    if (verbose) {
      cat("Updating seq, bias, and od.\n\n")
    }
    fout_seq <- foreach::foreach(index = seq_len(nsnps), .combine = cbind,
                                 .export = c("obj_for_eps",
                                             "grad_for_eps")) %dopar% {
      oout <- stats::optim(par       = c(seq[index], bias[index], od[index]),
                           fn        = obj_for_eps,
                           gr        = grad_for_eps,
                           method    = "L-BFGS-B",
                           lower     = rep(10^-6, 3),
                           upper     = c(1 - 10^-6, Inf, 1 - 10^-6),
                           control   = list(fnscale = -1, maxit = 10),
                           refvec    = refmat[, index],
                           sizevec   = sizemat[, index],
                           ploidy    = ploidy,
                           mean_bias = mean_bias,
                           var_bias  = var_bias,
                           mean_seq  = mean_seq,
                           var_seq   = var_seq,
                           mean_od   = 0,
                           var_od    = Inf,
                           wmat      = warray[, index, ])
      oout$par
      }
    seq  <- fout_seq[1, ]
    bias <- fout_seq[2, ]
    od   <- fout_seq[3, ]

    lbeta_array <- compute_all_log_bb(refmat = refmat,
                                      sizemat = sizemat,
                                      ploidy = ploidy,
                                      seq = seq,
                                      bias = bias,
                                      od = od)

    ## Calculate objective ---------------------------------------------------------------------------------------------------------
    obj <- elbo(warray      = warray,
                lbeta_array = lbeta_array,
                cor_inv     = cor_inv,
                postmean    = postmean,
                postvar     = postvar,
                bias        = bias,
                seq         = seq,
                mean_bias   = mean_bias,
                var_bias    = var_bias,
                mean_seq    = mean_seq,
                var_seq     = var_seq,
                ploidy      = ploidy)

    ## stopping criteria -----------------------------------------------------------------------------------------------------------
    err <- abs(obj_old / obj) - 1

    cat(" iteration:", iter, "\n",
        "objective:", obj, "\n",
        "      err:", err, "\n\n")

    iter <- iter + 1
  }

  ## Close cluster if used parallelization ----------------------------
  if (num_core > 1) {
    parallel::stopCluster(cl)
  }

  map_dosage <- apply(warray, c(1, 2), which.max) - 1
  maxpostprob <- apply(warray, c(1, 2), max)

  return_list               <- list()
  return_list$map_dosage    <- map_dosage
  return_list$maxpostprob   <- maxpostprob
  return_list$postprob      <- warray
  return_list$seq           <- seq
  return_list$bias          <- bias
  return_list$od            <- od
  return_list$allele_freq   <- allele_freq
  return_list$inbreeding    <- inbreeding
  return_list$cor_mat       <- cor_mat
  return_list$postmean      <- postmean
  return_list$postvar       <- postvar
  return_list$input$refmat  <- refmat
  return_list$input$sizemat <- sizemat
  return_list$input$ploidy  <- ploidy
  return_list$obj           <- obj

  class(return_list) <- "mupdog"
  return(return_list)
}

#' Update the underlying correlation matrix.
#'
#' @param postmean The matrix of posterior means. The rows index the
#'     individuals and the columns index the SNPs.
#' @param postvar The matrix of posterior variances. The rows index
#'     the individuals and the columns index the SNPs.
#'
#' @return A symmetric matrix of numerics. The update of the underlying
#'     correlation matrix.
#'
#' @author David Gerard
#'
#' @keywords internal
update_R <- function(postmean, postvar) {
  # stats::cov2cor(tcrossprod(postmean) + diag(rowSums(postvar)))
  (tcrossprod(postmean) + diag(rowSums(postvar))) / ncol(postmean)
}



#' Tests if its argument is a mupdog object.
#'
#' @param x Anything.
#'
#' @return A logical. \code{TRUE} if \code{x} is a mupdog object, and \code{FALSE} otherwise.
#'
#' @author David Gerard
#'
#' @examples
#' is.mupdog("anything")
#' # FALSE
#'
#' @export
#'
is.mupdog <- function(x) inherits(x, "mupdog")
