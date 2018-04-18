## Functions for ashdog and flexdog


#' Flexible genotyping for autopolyploids from next-generation sequencing data.
#'
#' This function will genotype polyploid individuals from next generation
#' sequencing (NGS) data while assuming the genotype distribution is either
#' unimodal (\code{model = "ash"}), generically any categorical
#' distribution (\code{model = "flex"}), binomial as a result of assuming
#' the population is in Hardy-Weinberg equlibrium (\code{model = "hw"}),
#' an overdispersed binomial (\code{model = "bb"}),
#' results from the individuals being siblings from either an F1 cross
#' (\code{model = "f1"}) or an S1 cross (\code{model = "s1"}), or is a discrete
#' uniform distribution (\code{model = "uniform"}).
#' It does this while accounting for
#' allele bias, overdispersion, and sequencing error.
#'
#' You might think a good default is \code{model = "uniform"} because it is
#' somehow an "uninformative prior." But it is very informative and tends to
#' work horribly in practice. I include it as an option only for completeness.
#'
#' You might also think that allowing the genotype distribution to be anything
#' would work well (\code{model = "flex"}). But this tends to overfit the data
#' and get stuck in local modes during optimization.
#'
#' A good default is either \code{model = "hw"}
#' if you have a small number of individuals (say, \eqn{<500}),
#' or \code{model = "ash"} if you have a lot of inidivuals (say, \eqn{>500}).
#' If the individuals are all siblings, then \code{model = "f1"} or
#' \code{model = "s1"} would work better. If the relatedness pattern between
#' individuals is more complicated, then I might recommend trying out
#' \code{\link{mupdog}}.
#'
#' Prior to using \code{flexdog}, during the read-mapping step,
#' you could try to get rid of allelic bias by
#' using WASP (\url{https://doi.org/10.1101/011221}). If you are successful
#' in removing the allelic bias (because its only source was the read-mapping
#' step), then you could set \code{update_bias = FALSE}. You can visually
#' inspect SNPs for bias by using \code{\link[updog]{plot_geno}} from the
#' \code{updog} package.
#'
#' @param refvec A vector of counts of reads with the reference allele.
#' @param sizevec A vector of total counts.
#' @param ploidy The ploidy of the species. Assumed to be the same for each
#'     individual.
#' @param model What form should the prior take? Should the genotype
#'     distribution be unimodal (\code{"ash"}), generically
#'     any categorical distribution (\code{"flex"}), binomial as a
#'     result of assuming Hardy-Weinberg equilibrium (\code{"hw"}),
#'     an overdispersed binomial (\code{"bb"}),
#'     a convolution of hypergeometics as a result that the population
#'     consists of either an F1 cross (\code{"f1"}) or an S1
#'     cross (\code{"s1"}), or fixed at a discrete uniform
#'     (\code{"uniform"})?
#' @param verbose Should we output more (\code{TRUE}) or less
#'     (\code{FALSE})?
#' @param mean_bias The prior mean of the log-bias.
#' @param var_bias The prior variance of the log-bias.
#' @param mean_seq The prior mean of the logit of the sequencing
#'     error rate.
#' @param var_seq The prior variance of the logit of the sequencing
#'     error rate.
#' @param seq The starting value of the sequencing error rate.
#' @param bias The starting value of the bias.
#' @param od The starting value of the overdispersion parameter.
#' @param mode The mode if \code{model = "ash"}. If not provided,
#'     \code{flexdog} will estimate the mode. This is the starting point of
#'     the allele frequency if \code{model = "hw"}. This should be
#'     \code{NULL} for all other options of \code{model}.
#' @param itermax The maximum number of EM iterations to run for each mode
#'     (if \code{model = "ash"}) or the total number of EM iterations to
#'     run (if \code{model = "flex"}, \code{model = "hw"},
#'     \code{model = "f1"}, \code{model = "s1"}, or \code{model = "uniform"}).
#' @param tol The tolerance stopping criterion. The EM algorithm will stop
#'     if the difference in the log-likelihoods between two consecutive
#'     iterations is less than \code{tol}.
#' @param use_cvxr A logical. If \code{model = "ash"}, then do you want
#'     to use the EM algorithm
#'     (\code{FALSE}) or a convex optimization program using
#'     the package CVXR \code{TRUE}?
#'     Only available if CVXR is installed.
#' @param update_bias A logical. Should we update \code{bias}
#'     (\code{TRUE}), or not (\code{FALSE})?
#' @param update_seq A logical. Should we update \code{seq}
#'     (\code{TRUE}), or not (\code{FALSE})?
#' @param update_od A logical. Should we update \code{od}
#'     (\code{TRUE}), or not (\code{FALSE})?
#' @param fs1_alpha Either \code{"optim"} or the value at which to fix
#'     the mixing proportion when \code{model = "f1"} or
#'     \code{model = "s1"}. If \code{optim}, then we optimize over
#'     the mixing proportion each iteration with possible values between
#'     \code{10^-8} and \code{10^-3}. If you fix it, I would recommend some small
#'     value such at \code{10^-3}.
#' @param p1ref The reference counts for the first parent if
#'     \code{model = "f1"}, or for the only parent if \code{model = "s1"}.
#' @param p1size The total counts for the first parent if
#'     \code{model = "f1"}, or for the only parent if \code{model = "s1"}.
#' @param p2ref The reference counts for the second parent if
#'     \code{model = "f1"}.
#' @param p2size The total counts for the second parent if
#'     \code{model = "f1"}.
#' @param ashpen The penalty to put on the unimodal prior.
#'
#' @return An object of class \code{flexdog}, which consists
#'     of a list with some or all of the following elements:
#' \describe{
#'   \item{\code{bias}}{The estimated bias parameter.}
#'   \item{\code{seq}}{The estimated sequencing error rate.}
#'   \item{\code{od}}{The estimated overdispersion parameter.}
#'   \item{\code{num_iter}}{The number of EM iterations ran. You should
#'       be wary if this equals \code{itermax}.}
#'   \item{\code{llike}}{The maximum marginal log-likelihood.}
#'   \item{\code{postmat}}{A matrix of posterior probabilities of each
#'       genotype for each individual. The rows index the individuals
#'       and the columns index the allele dosage.}
#'   \item{\code{gene_dist}}{The estimated genotype distribution. The
#'       \code{i}th element is the proportion of individuals with
#'       genotype \code{i-1}.}
#'   \item{\code{par}}{A list of the final estimates of the parameters
#'       of the genotype distribution. If \code{model = "hw"} then
#'       this will consist of \code{alpha}, the allele frequency.
#'       If \code{model = "f1"} or \code{model = "s1"} then this will
#'       consist of the parent genotype(s) and the mixing proportion
#'       with the discrete uniform (also called \code{alpha}).}
#'   \item{\code{geno}}{The posterior mode genotype.}
#'   \item{\code{maxpostprob}}{The maximum posterior probability.}
#'   \item{\code{postmean}}{The posterior mean genotype.}
#'   \item{\code{input$refvec}}{The value of \code{refvec} provided by
#'       the user.}
#'   \item{\code{input$sizevec}}{The value of \code{sizevec} provided
#'       by the user.}
#'   \item{\code{input$ploidy}}{The value of \code{ploidy} provided
#'       by the user.}
#'   \item{\code{input$model}}{The value of \code{model} provided by
#'       the user.}
#'   \item{\code{input$p1ref}}{The value of \code{p1ref} provided by the user.}
#'   \item{\code{input$p1size}}{The value of \code{p1size} provided by the user.}
#'   \item{\code{input$p2ref}}{The value of \code{p2ref} provided by the user.}
#'   \item{\code{input$p2size}}{The value of \code{p2size} provided by the user.}
#'   \item{\code{prop_mis}}{The posterior proportion of individuals
#'       misclassified.}
#' }
#'
#' @author David Gerard
#'
#' @export
flexdog <- function(refvec,
                    sizevec,
                    ploidy,
                    model       = c("hw", "bb", "ash", "f1", "s1", "flex", "uniform"),
                    verbose     = TRUE,
                    mean_bias   = 0,
                    var_bias    = 0.7 ^ 2,
                    mean_seq    = -4.7,
                    var_seq     = 1,
                    seq         = 0.005,
                    bias        = 1,
                    od          = 0.001,
                    update_bias = TRUE,
                    update_seq  = TRUE,
                    update_od   = TRUE,
                    mode        = NULL,
                    use_cvxr    = FALSE,
                    itermax     = 200,
                    tol         = 10 ^ -4,
                    fs1_alpha   = 10 ^ -3,
                    ashpen      = 0,
                    p1ref       = NULL,
                    p1size      = NULL,
                    p2ref       = NULL,
                    p2size      = NULL) {

  ## Check input -----------------------------------------------------
  model <- match.arg(model)
  if (model == "uniform") {
    warning("flexdog: Using model = 'uniform' is almost always a bad idea.\nTry model = 'hw' if you have data from a population study.")
  }

  assertthat::are_equal(length(refvec), length(sizevec))
  assertthat::assert_that(all(sizevec >= refvec, na.rm = TRUE))
  assertthat::assert_that(all(refvec >= 0, na.rm = TRUE))
  assertthat::are_equal(1, length(verbose), length(mean_bias),
                        length(var_bias), length(mean_seq),
                        length(var_seq), length(seq),
                        length(bias), length(od))
  assertthat::assert_that(var_bias > 0)
  assertthat::assert_that(var_seq > 0)
  assertthat::assert_that(seq >= 0, seq <= 1)
  assertthat::assert_that(bias > 0)
  assertthat::assert_that(od >= 0, od <= 1)
  assertthat::assert_that(is.logical(verbose))
  assertthat::are_equal(ploidy %% 1, 0)
  assertthat::assert_that(ploidy > 0)
  assertthat::assert_that(tol > 0)
  assertthat::are_equal(itermax %% 1, 0)
  assertthat::assert_that(itermax > 0)
  assertthat::assert_that(is.logical(use_cvxr))
  assertthat::assert_that(is.logical(update_bias))
  assertthat::assert_that(is.logical(update_seq))
  assertthat::assert_that(is.logical(update_od))
  assertthat::assert_that(ashpen >= 0)

  ## check fs1_alpha -----------------------------------------------
  if (is.numeric(fs1_alpha)) {
    assertthat::are_equal(length(fs1_alpha), 1)
    assertthat::assert_that(fs1_alpha <= 1, fs1_alpha >= 0)
  } else {
    if (fs1_alpha != "optim") {
      stop("flexdog: fs1_alpha either needs to be 'optim' or a numeric between 0 and 1.")
    }
  }

  ## Check p1ref, p2ref, p1size, p2size ----------------------------
  if ((!is.null(p1ref) | !is.null(p1size)) & (model != "f1" & model != "s1")) {
    stop("flexdog: if model is not 'f1' or 's1', then p1ref and p1size both need to be NULL.")
  }
  if ((!is.null(p2ref) | !is.null(p2size)) & (model != "f1")) {
    stop("flexdog: if model is not 'f1', then p2ref and p2size both need to be NULL.")
  }
  if ((is.null(p1ref) & !is.null(p1size)) | (!is.null(p1ref) & is.null(p1size))) {
    stop("flexdog: p1ref and p1size either need to be both NULL or both non-NULL.")
  }
  if ((is.null(p2ref) & !is.null(p2size)) | (!is.null(p2ref) & is.null(p2size))) {
    stop("flexdog: p1ref and p1size either need to be both NULL or both non-NULL.")
  }
  if (!is.null(p1ref)) {
    assertthat::assert_that(p1ref >= 0, p1size >= p1ref)
  }
  if (!is.null(p2ref)) {
    assertthat::assert_that(p2ref >= 0, p2size >= p2ref)
  }

  ## check and set mode under various models -----------------------
  if (!is.null(mode) & model == "flex") {
    stop('flexdog: `model` cannot equal `"flex"` when `mode` is specified.')
  } else if (is.null(mode) & model == "flex") {
    mode_vec <- 0
  } else if (!is.null(mode) & (model == "f1" | model == "s1" | model == "uniform")) {
    stop('flexdog: `model` cannot equal `"f1" or "s1"` when `mode` is specified.')
  } else if (is.null(mode) & (model == "f1" | model == "s1" | model == "uniform")) {
    mode_vec <- mean(refvec / sizevec, na.rm = TRUE) ## just to initialize pivec
  } else if (!is.null(mode) & model == "ash") {
    assertthat::are_equal(length(mode), 1)
    mode_vec <- mode
  } else if (is.null(mode) & model == "ash") {
    mode_vec <- (0:(ploidy - 1)) + 0.5
  } else if (is.null(mode) & (model == "hw" | model == "bb")) {
    mode_vec <- mean(refvec / sizevec, na.rm = TRUE)
  } else if (!is.null(mode) & (model == "hw" | model == "bb")) {
    if (any((mode < 0) | (mode > 1))) {
      stop('If model = "hw" or model = "bb" then `mode` should be between 0 and 1.\nIt is the initialization of the allele frequency.')
    }
  } else {
    stop("flexdog: Checking mode. How did you get here?")
  }

  ## Deal with missingness in sizevec and refvec -----------------------
  not_na_vec  <- !(is.na(refvec) | is.na(sizevec))
  refvec      <- refvec[not_na_vec]
  sizevec     <- sizevec[not_na_vec]

  ## Some variables needed to run EM ---------------------------
  control <- list() ## will contain parameters used to update pivec
  boundary_tol <- 10 ^ -6 ## smallest values of seq, bias, od,
                          ## how close can od and seq get to 1.
  if (od < boundary_tol) {
    od <- boundary_tol
  }
  if (seq < boundary_tol) {
    seq <- boundary_tol
  }
  if (bias < boundary_tol) {
    bias <- boundary_tol
  }
  if (od > 1 - boundary_tol) {
    od <- 1 - boundary_tol
  }
  if (seq > 1 - boundary_tol) {
    seq <- 1 - boundary_tol
  }

  ## Run EM for each mode in `mode_vec` -----------------------
  return_list <- list(llike = -Inf)
  for (em_index in 1:length(mode_vec)) {
    mode <- mode_vec[em_index]

    if (verbose) {
      cat("Mode:", mode, "\n")
    }

    ## Get inner weight vec only once
    ## Used in convex optimization program
    ## assign use_cvxr
    if (model == "ash") {
      control$inner_weights <- get_inner_weights(ploidy = ploidy, mode = mode)
      control$use_cvxr      <- use_cvxr
      control$lambda        <- ashpen
    } else if (model == "f1" | model == "s1") {
      control$qarray    <- updog::get_q_array(ploidy = ploidy)
      control$fs1_alpha <- fs1_alpha
    } else if (model == "bb") {
      control$alpha <- mode ## initialize allele frequency for bb
      control$tau   <- boundary_tol ## initialize od for bb
    }

    ## Initialize pivec so that two modes have equal prob if model = "ash".
    ##     Uniform if model = "flex".
    pivec <- initialize_pivec(ploidy = ploidy, mode = mode, model = model)
    assertthat::are_equal(sum(pivec), 1)
    control$pivec <- pivec ## initialization for unimodal optimization

    probk_vec <- get_probk_vec(pivec = pivec, model = model, mode = mode)
    assertthat::are_equal(sum(probk_vec), 1)

    ## Run EM ----------------------------------------
    iter_index  <- 1
    err         <- tol + 1
    llike       <- -Inf
    while (err > tol & iter_index <= itermax) {
      llike_old <- llike

      ## E-step ----------------------
      wik_mat <- get_wik_mat(probk_vec = probk_vec, refvec = refvec,
                             sizevec = sizevec, ploidy = ploidy,
                             seq = seq, bias = bias, od = od)

      ## Update seq, bias, and od ----
      oout <- stats::optim(par         = c(seq, bias, od),
                           fn          = obj_for_eps,
                           gr          = grad_for_eps,
                           method      = "L-BFGS-B",
                           lower       = rep(boundary_tol, 3),
                           upper       = c(1 - boundary_tol, Inf,
                                           1 - boundary_tol),
                           control     = list(fnscale = -1, maxit = 20),
                           refvec      = refvec,
                           sizevec     = sizevec,
                           ploidy      = ploidy,
                           mean_bias   = mean_bias,
                           var_bias    = var_bias,
                           mean_seq    = mean_seq,
                           var_seq     = var_seq,
                           wmat        = wik_mat,
                           update_seq  = update_seq,
                           update_bias = update_bias,
                           update_od   = update_od)
      seq  <- oout$par[1]
      bias <- oout$par[2]
      od   <- oout$par[3]

      ## if F1 or S1, update betabinomial log-likelihood of parent counts
      if (model == "f1") {
        if (!is.null(p1ref)) {
          xi_vec <- xi_fun(p = 0:ploidy / ploidy, eps = seq, h = bias)
          control$p1_lbb <- dbetabinom(x    = rep(p1ref, ploidy + 1),
                                       size = rep(p1size, ploidy + 1),
                                       mu   = xi_vec,
                                       rho  = od,
                                       log  = TRUE)
        }
        if (!is.null(p2ref)) {
          xi_vec <- xi_fun(p = 0:ploidy / ploidy, eps = seq, h = bias)
          control$p2_lbb <- dbetabinom(x    = rep(p2ref, ploidy + 1),
                                       size = rep(p2size, ploidy + 1),
                                       mu   = xi_vec,
                                       rho  = od,
                                       log  = TRUE)
        }
      } else if (model == "s1") {
        if (!is.null(p1ref)) {
          xi_vec <- xi_fun(p = 0:ploidy / ploidy, eps = seq, h = bias)
          control$p1_lbb <- dbetabinom(x    = rep(p1ref, ploidy + 1),
                                       size = rep(p1size, ploidy + 1),
                                       mu   = xi_vec,
                                       rho  = od,
                                       log  = TRUE)
        }
      } else {
        ## do nothing
      }

      ## Update pivec ----------------
      weight_vec <- colSums(wik_mat)
      fupdate_out <- flex_update_pivec(weight_vec = weight_vec,
                                       model      = model,
                                       control    = control)
      pivec <- fupdate_out$pivec
      control$pivec <- pivec ## initial condition for unimodal optimization

      if (model == "bb") { ## update alpha and tau in control
        control$alpha <- fupdate_out$par$alpha
        control$tau   <- fupdate_out$par$tau
      }


      ## Update probk_vec -----------------------------------------------
      pivec[pivec < 0] <- 0
      pivec[pivec > 1] <- 1
      probk_vec <- get_probk_vec(pivec = pivec, model = model, mode = mode)

      ## Calculate likelihood and update stopping criteria --------------
      llike <- flexdog_obj(probk_vec = probk_vec,
                           refvec    = refvec,
                           sizevec   = sizevec,
                           ploidy    = ploidy,
                           seq       = seq,
                           bias      = bias,
                           od        = od,
                           mean_bias = mean_bias,
                           var_bias  = var_bias,
                           mean_seq  = mean_seq,
                           var_seq   = var_seq)

      if (model == "ash" & !use_cvxr) { ## add small penalty if "ash"
        llike <- llike + ashpen_fun(lambda = ashpen, pivec = pivec)
      }

      err        <- abs(llike - llike_old)
      iter_index <- iter_index + 1

      if (llike < llike_old - 10 ^ -5) {
        warning(paste0("flexdog: likelihood not increasing.\nDifference is: ",
                       llike - llike_old))
        if (verbose) {
          cat("\nindex: ", iter_index, "\n")
          cat("llike: ", llike, "\n")
          cat("pivec: ", paste0("c(", paste0(pivec, collapse = ", "), ")"), "\n")
          cat("weight_vec: ",  paste0("c(", paste0(weight_vec, collapse = ", "), ")"), "\n\n")
        }
      }
    }

    if (verbose) {
      cat("llike:", llike, "\n\n")
    }

    ## Check which mode has the highest likelihood --------------------
    temp_list <- list(bias      = bias,
                      seq       = seq,
                      od        = od,
                      num_iter  = iter_index,
                      llike     = llike,
                      postmat   = wik_mat,
                      gene_dist = probk_vec,
                      par       = fupdate_out$par)
    if (temp_list$llike > return_list$llike) {
      return_list <- temp_list
    }
  }

  ## Summaries --------------------------------------------------------
  return_list$geno          <- apply(return_list$postmat, 1, which.max) - 1
  return_list$maxpostprob   <- return_list$postmat[cbind(1:nrow(return_list$postmat), return_list$geno + 1)]
  return_list$postmean      <- c(return_list$postmat %*% 0:ploidy)
  return_list$input$refvec  <- refvec
  return_list$input$sizevec <- sizevec
  return_list$input$ploidy  <- ploidy
  return_list$input$model   <- model
  return_list$input$p1ref   <- p1ref
  return_list$input$p1size  <- p1size
  return_list$input$p2ref   <- p2ref
  return_list$input$p2size  <- p2size
  return_list$prop_mis      <- 1 - mean(return_list$maxpostprob)


  ## Add back missingness ---------------------------------------------
  temp                     <- rep(NA, length = length(not_na_vec))
  temp[not_na_vec]         <- refvec
  return_list$input$refvec <- temp

  temp                      <- rep(NA, length = length(not_na_vec))
  temp[not_na_vec]          <- sizevec
  return_list$input$sizevec <- temp

  temp             <- rep(NA, nrow = length(not_na_vec))
  temp[not_na_vec] <- return_list$geno
  return_list$geno <- temp

  temp                    <- rep(NA, nrow = length(not_na_vec))
  temp[not_na_vec]        <- return_list$maxpostprob
  return_list$maxpostprob <- temp

  temp                 <- rep(NA, nrow = length(not_na_vec))
  temp[not_na_vec]     <- return_list$postmean
  return_list$postmean <- temp

  temp                <- matrix(NA, nrow = length(not_na_vec),
                                ncol = ncol(return_list$postmat))
  temp[not_na_vec, ]  <- return_list$postmat
  return_list$postmat <- temp

  ## Set class to flexdog ---------------------------------------------
  class(return_list) <- "flexdog"

  return(return_list)
}

#' Draw a genotype plot from the output of \code{\link{flexdog}}.
#'
#' @inherit plot.mupdog description details
#'
#' @param x A \code{flexdog} object.
#' @param ... Not used.
#'
#' @author David Gerard
#'
#' @seealso
#' \describe{
#'   \item{\code{\link[updog]{plot_geno}}}{The underlying plotting function.}
#'   \item{\code{\link{flexdog}}}{Creates a \code{flexpdog} object.}
#' }
#'
#' @export
plot.flexdog <- function(x, ...) {
  assertthat::assert_that(is.flexdog(x))
  pl <- updog::plot_geno(ocounts   = x$input$refvec,
                         osize     = x$input$sizevec,
                         ploidy    = x$input$ploidy,
                         ogeno     = x$geno,
                         seq_error = x$seq,
                         bias_val  = x$bias,
                         prob_ok   = x$maxpostprob) +
    ggplot2::guides(alpha=ggplot2::guide_legend(title = "maxpostprob"))
  return(pl)
}

#' Tests if an argument is a \code{flexdog} object.
#'
#' @param x Anything.
#'
#' @return A logical. \code{TRUE} if \code{x} is a \code{flexdog} object, and \code{FALSE} otherwise.
#'
#' @author David Gerard
#'
#' @export
is.flexdog <- function(x) {
  inherits(x, "flexdog")
}


#' Initialize \code{pivec} for \code{\link{flexdog}} EM algorithm.
#'
#' The key idea here is choosing the pi's so that the two modes
#' have equal probability.
#'
#' @inheritParams flexdog
#'
#' @seealso \code{\link{flexdog}} for where this is used.
#'
#' @author David Gerard
initialize_pivec <- function(ploidy, mode, model = c("hw", "bb", "ash", "f1", "s1", "flex", "uniform")) {
  assertthat::are_equal(1, length(ploidy), length(mode))
  assertthat::are_equal(ploidy %% 1, 0)

  model <- match.arg(model)
  if (model == "flex" | model == "uniform") {
    pivec <- rep(x = 1 / (ploidy + 1), length = ploidy + 1)
  } else if (model == "ash") {
    init_type <- "equi"
    if (init_type == "equi") { ## equimodal
      floor_mode <- floor(mode)
      ceil_mode  <- ceiling(mode)
      d <- sum(1 / (ceil_mode:ploidy - ceil_mode + 1)) /
        (sum(1 / (floor_mode - 0:floor_mode + 1)))
      second_half <- 1 / (ploidy - ceil_mode + 1 + (floor_mode + 1) * d)
      first_half <- second_half * d
      ## assertthat::are_equal(first_half, (1 - (ploidy - ceil_mode + 1) * second_half) / (floor_mode + 1))
      pivec <- c(rep(first_half, length = floor_mode + 1),
                 rep(second_half, length = ploidy - ceil_mode + 1))
    } else if (init_type == "bin") { ## binomial
      pvec_init <- stats::dbinom(x = 0:ploidy, size = ploidy,
                                 prob = floor(mode) / ploidy)
      pivec <- get_uni_rep(pvec_init)$pivec + 10 ^-6
      pivec <- pivec / sum(pivec)
    }
  } else if (model == "hw" | model == "f1" | model == "s1" | model == "bb") {
    if (mode < 0 | mode > 1) {
      stop('initialize_pivec: when model = "hw", mode should be between 0 and 1.\n It is the initialization of the allele frequency.')
    }
    pivec <- stats::dbinom(x = 0:ploidy, size = ploidy, prob = mode)
  } else {
    stop("initialize_pivec: How did you get here?")
  }
  return(pivec)
}


#' Update the distribution of genotypes from various models.
#'
#' @param weight_vec \code{colSums(wik_mat)} from \code{\link{flexdog}}.
#'     This is the sum of current posterior probabilities of each individual
#'     having genotype k.
#' @param model What model are we assuming.
#' @param control A list of anything else needed to be passed.
#'     E.g. if \code{model = "ash"},
#'     then \code{inner_weights} needs to be passed through \code{control}
#'     (see \code{\link{get_inner_weights}} for how to get this matrix).
#'
#' @return A list with the following elements
#' \describe{
#'   \item{\code{pivec}}{The estimate of the genotype distribution.}
#'   \item{\code{par}}{A list of estimated parameters. An empty list if the model does not contain any parameters other than \code{pivec}.}
#' }
#'
#' @author David Gerard
flex_update_pivec <- function(weight_vec, model = c("hw", "bb", "ash", "f1", "s1", "flex", "uniform"), control) {
  ## Check input -------------------------------
  ploidy <- length(weight_vec) - 1
  model <- match.arg(model)
  assertthat::are_equal(sum(weight_vec), 1)
  assertthat::assert_that(is.list(control))

  ## Get pivec ---------------------------------
  return_list <- list()
  if (model == "flex") {
    pivec <- weight_vec / sum(weight_vec)
    return_list$pivec <- pivec
    return_list$par <- list()
  } else if (model == "ash") {
    if (is.null(control$inner_weights)) {
      stop('flex_update_pivec: control$inner_weights cannot be NULL when model = "ash"')
    }
    if (!control$use_cvxr) {
      pivec <- uni_em(weight_vec = weight_vec,
                      lmat       = control$inner_weights,
                      pi_init    = control$pivec,
                      itermax    = 200,
                      obj_tol    = 10 ^ -4,
                      lambda     = control$lambda)
    } else if (control$use_cvxr & requireNamespace("CVXR", quietly = TRUE)) {
      cv_pi <- CVXR::Variable(1, ploidy + 1)
      obj   <- sum(t(weight_vec) * log(cv_pi %*% control$inner_weights))
      prob  <- CVXR::Problem(CVXR::Maximize(obj),
                             constraints = list(sum(cv_pi) == 1,
                                                cv_pi >= 0))
      result <- solve(prob)
      result$value
      pivec <- c(result$getValue(cv_pi))
    } else {
      stop("flex_update_pivec: CVXR not installed but use_cvxr = TRUE.")
    }
    return_list$pivec <- pivec
    return_list$par <- list()
  } else if (model == "hw") {
    alpha <- sum(0:ploidy * weight_vec) / (ploidy * sum(weight_vec))
    pivec <- stats::dbinom(x = 0:ploidy, size = ploidy, prob = alpha)
    return_list$pivec     <- pivec
    return_list$par       <- list()
    return_list$par$alpha <- alpha
  } else if (model == "f1") {
    optim_best       <- list()
    optim_best$value <- -Inf
    for (i in 0:ploidy) { ## parent 1
      for (j in 0:ploidy) { ## parent 2
        pvec <- control$qarray[i + 1, j + 1, ]
        if (control$fs1_alpha == "optim") {
          optim_out <- stats::optim(par = 0.01,
                                    fn = f1_obj,
                                    method = "Brent",
                                    control = list(fnscale = -1),
                                    upper = 0.001,
                                    lower = 10 ^ -8,
                                    pvec = pvec,
                                    weight_vec = weight_vec)
        } else {
          optim_out <- list()
          optim_out$value <- f1_obj(alpha = control$fs1_alpha,
                                    pvec = pvec,
                                    weight_vec = weight_vec)
          optim_out$par   <- control$fs1_alpha
        }
        if (!is.null(control$p1_lbb)) {
          optim_out$value <- optim_out$value + control$p1_lbb[i + 1]
        }
        if (!is.null(control$p2_lbb)) {
          optim_out$value <- optim_out$value + control$p2_lbb[j + 1]
        }
        if (optim_out$value > optim_best$value) {
          optim_best <- optim_out
          optim_best$ell1 <- i
          optim_best$ell2 <- j
        }
      }
    }
    return_list$pivec <- (1 - optim_best$par) * control$qarray[optim_best$ell1 + 1, optim_best$ell2 + 1, ] +
      optim_best$par / (ploidy + 1)
    return_list$par <- list()
    return_list$par$p1geno <- optim_best$ell1
    return_list$par$p2geno <- optim_best$ell2
    return_list$par$alpha  <- optim_best$par
  } else if (model == "s1") {
    optim_best       <- list()
    optim_best$value <- -Inf
    for (i in 0:ploidy) { ## parent
      pvec <- control$qarray[i + 1, i + 1, ]
      if (control$fs1_alpha == "optim") {
        optim_out <- stats::optim(par = 0.01,
                                  fn = f1_obj,
                                  method = "Brent",
                                  control = list(fnscale = -1),
                                  upper = 0.001,
                                  lower = 10 ^ -8,
                                  pvec = pvec,
                                  weight_vec = weight_vec)
      } else {
        optim_out <- list()
        optim_out$value <- f1_obj(alpha = control$fs1_alpha,
                                  pvec = pvec,
                                  weight_vec = weight_vec)
        optim_out$par   <- control$fs1_alpha
      }
      if (!is.null(control$p1_lbb)) {
        optim_out$value <- optim_out$value + control$p1_lbb[i + 1]
      }
      if (optim_out$value > optim_best$value) {
        optim_best <- optim_out
        optim_best$ell <- i
      }
    }
    return_list$pivec <- (1 - optim_best$par) * control$qarray[optim_best$ell + 1, optim_best$ell + 1, ] +
      optim_best$par / (ploidy + 1)
    return_list$par <- list()
    return_list$par$pgeno <- optim_best$ell
    return_list$par$alpha  <- optim_best$par
  } else if (model == "uniform") {
    return_list$pivec <- rep(x = 1 / (ploidy + 1), length = ploidy + 1)
    return_list$par <- list()
  } else if (model == "bb") {
    boundary_val <- 10 ^ -8
    optim_out <- stats::optim(par        = c(control$alpha, control$tau),
                              fn         = obj_for_weighted_lbb,
                              gr         = grad_for_weighted_lbb,
                              method     = "L-BFGS-B",
                              lower      = c(boundary_val, boundary_val),
                              upper      = c(1 - boundary_val, 1 - boundary_val),
                              weight_vec = weight_vec,
                              ploidy     = ploidy,
                              control = list(fnscale = -1))

    return_list$pivec <- dbetabinom(x    = 0:ploidy,
                                    size = ploidy,
                                    mu   = optim_out$par[1],
                                    rho  = optim_out$par[2],
                                    log  = FALSE)
    return_list$par <- list()
    return_list$par$alpha <- optim_out$par[1]
    return_list$par$tau   <- optim_out$par[2]
  } else {
    stop("flex_update_pivec: how did you get here?")
  }
  return(return_list)
}



#' Get the representation of a discrete unimodal probability distribution.
#'
#' NB: In \code{get_uni_rep}, we count the mode starting at 1.
#'     In \code{\link{get_probk_vec}}, we count the mode starting at 0.
#'
#' @param probvec A probability vector. It assumes the probabilities
#'     are ordered according to the ordering of the discrete set.
#'
#' @seealso \code{\link{get_probk_vec}} with option \code{model = "ash"}
#'     for the inverse of this function.
#'
#' @author David Gerard
#'
get_uni_rep <- function(probvec) {
  assertthat::are_equal(sum(probvec), 1)
  assertthat::assert_that(all(probvec >= 0))

  mode <- which.max(probvec)
  n    <- length(probvec)

  if (mode < n) {
    fak_vec <- c(mode:1, 1:(n - mode))
  } else {
    fak_vec <- n:1
  }

  if (mode < n) {
    aug_probvec <- c(0, probvec, 0)
    diffvec <- c(aug_probvec[2:(mode + 1)] - aug_probvec[1:mode],
                 aug_probvec[(mode + 2):(n + 1)] - aug_probvec[(mode + 3):(n + 2)])
  } else {
    aug_probvec <- c(0, probvec)
    diffvec <- aug_probvec[2:(n + 1)] - aug_probvec[1:n]
  }

  return_list <- list()
  return_list$pivec <- diffvec * fak_vec
  return_list$mode = mode + 0.5
  return(return_list)
}



#' Penalty on pivec used when \code{model = "ash"} in \code{\link{flexdog}}.
#'
#' @param lambda The penalty.
#' @param pivec The vector of mixing proportions for the component
#'     discrete uniform distributions.
#'
#' @author David Gerard
#'
ashpen_fun <- function(lambda, pivec) {
  lambda * sum(log(pivec))
}


