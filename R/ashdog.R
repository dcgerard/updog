## Functions for ashdog and flexdog


#' Flexible updog
#'
#' This function will genotype polyploid individuals from next generation
#' sequencing (NGS) data while assuming the genotype distribution is either
#' unimodal (\code{model = "ash"}), generically any categorical
#' distribution (\code{model = "flex"}), or Binomial as a result of assuming
#' the population is in Hardy-Weinberg equlibrium (\code{model = "hw"}).
#' It does this while accounting for
#' allele bias, overdispersion, and sequencing error.
#'
#' @param refvec A vector of counts of reads with the reference allele.
#' @param sizevec A vector of total counts.
#' @param ploidy The ploidy of the species.
#' @param model What form should the prior take? Should the genotype
#'     distribution be unimodal (\code{"ash"}), generically
#'     any categorical distribution (\code{"flex"}), or binomial as a
#'     result of assuming Hardy-Weinberg equilibrium (\code{"hw"})?
#' @param verbose Should we output more (\code{TRUE}) or less
#'     (\code{FALSE})?
#' @param mean_bias The prior mean of the log-bias.
#' @param var_bias The prior variance of the log-bias.
#' @param mean_seq The prior mean of the logit of the sequencing error rate.
#' @param var_seq The prior variance of the logit of the sequencing
#'     error rate.
#' @param seq The starting value of the sequencing error rate.
#' @param bias The starting value of the bias.
#' @param od The starting value of the overdispersion parameter.
#' @param mode The mode if \code{model = "ash"}. If not provided,
#'     \code{flexdog} will estimate the mode. This is the starting point of
#'     the allele frequency if \code{model = "hw"}.
#' @param itermax The maximum number of EM iterations to run for each mode
#'     (if \code{model = "ash"}) or the total number of EM iterations to
#'     run (if \code{model = "flex"} or \code{model = "hw"}).
#' @param tol The tolerance stopping criterion. The EM algorithm will stop
#'     if the difference in the log-likelihoods between two consecutive
#'     iterations is less than \code{tol}.
#'
#' @return An object of class \code{flexdog}, which consists
#'     of a list with some or all of the following elements:
#' \describe{
#'   \item{\code{bias}}{The estimated bias parameter.}
#'   \item{\code{seq}}{The estimated sequencing error rate.}
#'   \item{\code{od}}{The estimated overdispersion parameter.}
#'   \item{\code{num_iter}}{The number of EM iterations ran. You should be wary if this equals \code{itermax}.}
#'   \item{\code{llike}}{The maximum marginal log-likelihood.}
#'   \item{\code{postmat}}{A matrix of posterior probabilities of each genotype for each individual. The rows index the individuals and the columns index the allele dosage.}
#'   \item{\code{gene_dist}}{The estimated genotype distribution. The \code{i}th element is the proportion of individuals with genotype \code{i-1}.}
#'   \item{\code{par}}{A list of the final estimates of the parameters of the genotype distribution.}
#'   \item{\code{geno}}{The posterior mode genotype.}
#'   \item{\code{maxpostprob}}{The maximum posterior probability.}
#'   \item{\code{postmean}}{The posterior mean genotype.}
#'   \item{\code{input$refvec}}{The value of \code{refvec} provided by the user.}
#'   \item{\code{input$sizevec}}{The value of \code{sizevec} provided by the user.}
#'   \item{\code{input$ploidy}}{The value of \code{ploidy} provided by the user.}
#'   \item{\code{input$model}}{The value of \code{model} provided by the user.}
#'   \item{\code{prop_mis}}{The posterior proportion of individuals misclassified.}
#' }
#'
#' @author David Gerard
#'
#' @export
flexdog <- function(refvec,
                    sizevec,
                    ploidy,
                    model     = c("ash", "flex", "hw"),
                    verbose   = TRUE,
                    mean_bias = 0,
                    var_bias  = 0.7 ^ 2,
                    mean_seq  = -4.7,
                    var_seq   = 1,
                    seq       = 0.005,
                    bias      = 1,
                    od        = 0.001,
                    mode      = NULL,
                    itermax   = 200,
                    tol       = 10^-4) {

  ## Check input -----------------------------------------------------
  model <- match.arg(model)
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

  if (!is.null(mode) & model == "flex") {
    stop('flexdog: `model` cannot equal `"flex"` when `mode` is specified.')
  } else if (is.null(mode) & model == "flex") {
    mode_vec <- 0
  } else if (!is.null(mode) & model == "ash") {
    assertthat::are_equal(length(mode), 1)
    mode_vec <- mode
  } else if (is.null(mode) & model == "ash") {
    mode_vec <- (0:(ploidy - 1)) + 0.5
  } else if (is.null(mode) & model == "hw") {
    mode_vec <- mean(refvec / sizevec, na.rm = TRUE)
  } else if (!is.null(mode) & model == "hw") {
    if (any((mode < 0) | (mode > 1))) {
      stop('If model = "hw" then `mode` should be between 0 and 1.\nIt is the initialization of the allele frequency.')
    }
  } else {
    stop("flexdog: Checking mode. How did you get here?")
  }

  ## Deal with missingness in sizevec and refvec -----------------------
  not_na_vec  <- !(is.na(refvec) | is.na(sizevec))
  refvec      <- refvec[not_na_vec]
  sizevec     <- sizevec[not_na_vec]

  ## Some variables needed to run EM ---------------------------
  boundary_tol <- 10 ^ -6
  control      <- list()

  ## Run EM for each mode in `mode_vec` -----------------------
  return_list <- list(llike = -Inf)
  for (em_index in 1:length(mode_vec)) {
    mode <- mode_vec[em_index]

    if (verbose) {
      cat("Mode:", mode, "\n")
    }

    ## Get inner weight vec only once
    ## Used in convex optimization program
    if (model == "ash") {
      control$inner_weights <- get_inner_weights(ploidy = ploidy, mode = mode)
    }


    ## Initialize pivec so that two modes have equal prob if model = "ash".
    ##     Uniform if model = "flex".
    pivec <- initialize_pivec(ploidy = ploidy, mode = mode, model = model)
    assertthat::are_equal(sum(pivec), 1)

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
      oout <- stats::optim(par       = c(seq, bias, od),
                           fn        = obj_for_eps,
                           gr        = grad_for_eps,
                           method    = "L-BFGS-B",
                           lower     = rep(boundary_tol, 3),
                           upper     = c(1 - boundary_tol, Inf,
                                         1 - boundary_tol),
                           control   = list(fnscale = -1, maxit = 20),
                           refvec    = refvec,
                           sizevec   = sizevec,
                           ploidy    = ploidy,
                           mean_bias = mean_bias,
                           var_bias  = var_bias,
                           mean_seq  = mean_seq,
                           var_seq   = var_seq,
                           wmat      = wik_mat)
      seq  <- oout$par[1]
      bias <- oout$par[2]
      od   <- oout$par[3]

      ## Update pivec ----------------
      weight_vec <- colSums(wik_mat)
      fupdate_out <- flex_update_pivec(weight_vec = weight_vec,
                                       model      = model,
                                       control    = control)
      pivec <- fupdate_out$pivec

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
      err        <- abs(llike - llike_old)
      iter_index <- iter_index + 1

      if (llike < llike_old - 10 ^ -8) {
        warning(paste0("flexdog: likelihood not increasing.\nDifference is",
                       llike - llike_old))
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
    ggplot2::guides(alpha=ggplot2::guide_legend(title="maxpostprob"))
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
initialize_pivec <- function(ploidy, mode, model = c("ash", "flex", "hw")) {
  assertthat::are_equal(1, length(ploidy), length(mode))
  assertthat::are_equal(ploidy %% 1, 0)

  model <- match.arg(model)
  if (model == "flex") {
    pivec <- rep(x = 1 / (ploidy + 1), length = ploidy + 1)
  } else if ((mode <= 0) | (mode >= ploidy) | (mode %% 1 == 0)) {
    pivec <- rep(x = 1 / (ploidy + 1), length = ploidy + 1)
  } else if (model == "ash") {
    floor_mode <- floor(mode)
    ceil_mode  <- ceiling(mode)

    d <- sum(1 / (ceil_mode:ploidy - ceil_mode + 1)) /
      (sum(1 / (floor_mode - 0:floor_mode + 1)))

    second_half <- 1 / (ploidy - ceil_mode + 1 + (floor_mode + 1) * d)

    first_half <- second_half * d

    ## assertthat::are_equal(first_half, (1 - (ploidy - ceil_mode + 1) * second_half) / (floor_mode + 1))

    pivec <- c(rep(first_half, length = floor_mode + 1),
               rep(second_half, length = ploidy - ceil_mode + 1))
  } else if (model == "hw") {
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
flex_update_pivec <- function(weight_vec, model = c("ash", "flex", "hw"), control) {
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
    cv_pi <- CVXR::Variable(1, ploidy + 1)
    obj   <- sum(t(weight_vec) * log(cv_pi %*% control$inner_weights))
    prob  <- CVXR::Problem(CVXR::Maximize(obj),
                           constraints = list(sum(cv_pi) == 1,
                                              cv_pi >= 0))
    result <- solve(prob)
    result$value
    pivec <- c(result$getValue(cv_pi))
    return_list$pivec <- pivec
    return_list$par <- list()
  } else if (model == "hw") {
    alpha <- sum(0:ploidy * weight_vec) / (ploidy * sum(weight_vec))
    pivec <- stats::dbinom(x = 0:ploidy, size = ploidy, prob = alpha)
    return_list$pivec     <- pivec
    return_list$par       <- list()
    return_list$par$alpha <- alpha
  } else {
    stop("flex_update_pivec: how did you get here?")
  }
  return(return_list)
}


