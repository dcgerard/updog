## Functions for ashdog and flexdog


#' Flexible updog
#'
#' @param refvec A vector of counts of reads with reference allele.
#' @param sizevec A vector of total counts.
#' @param ploidy The ploidy of the species.
#' @param model What form should the prior take? Should the genotype distribution be unimodal
#'     (\code{"ash"}) or should be generically any categorical distribution (\code{"flex"})?
#' @param verbose Should we output more (\code{TRUE}) or less (\code{FALSE})?
#' @param mean_bias The prior mean of the log-bias.
#' @param var_bias The prior variance of the log-bias.
#' @param mean_seq The prior mean of the logit of the sequencing error rate.
#' @param var_seq The prior variance of the logit of the sequencing error rate.
#' @param seq The starting value of the sequencing error rate.
#' @param bias The starting value of the bias.
#' @param od The starting value of the overdispersion parameter.
#' @param mode The mode if \code{model = "ash"}. If not provided, \code{flexdog} will estimate the mode.
#'
#' @author David Gerard
#'
flexdog <- function(refvec,
                    sizevec,
                    ploidy,
                    model     = c("ash", "flex"),
                    verbose   = TRUE,
                    mean_bias = 0,
                    var_bias  = 1,
                    mean_seq  = -4.7,
                    var_seq   = 1,
                    seq       = 0.005,
                    bias      = 1,
                    od        = 0.001,
                    mode      = NULL) {

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

  if (!is.null(mode) & model == "flex") {
    stop('flexdog: `model` cannot equal `"flex"` when `mode` is specified.')
  } else if (is.null(mode) & model == "flex") {
    mode_vec <- -1
  } else if (!is.null(mode) & model == "ash") {
    assertthat::are_equal(length(mode), 1)
    assertthat::are_equal(mode %% 1, 0)
    mode_vec <- mode + 0.5
  } else if (is.null(mode) & model == "ash") {
    mode_vec <- (0:(ploidy - 1)) + 0.5
  }

  ## Run EM for each mode in `mode_vec`.
  for (em_index in 1:length(mode_vec)) {
    mode <- mode_vec[em_index]

    ## Initialize pivec so that two modes have equal prob if model = "ash". Uniform if model = "flex".
    pivec <- initialize_pivec(ploidy = ploidy, mode = mode, model = model)
    assertthat::are_equal(sum(pivec), 1)

    probk_vec <- get_probk_vec(pivec = pivec, model = model, mode = mode)
    assertthat::are_equal(sum(probk_vec), 1)
  }

}


#' Initialize \code{pivec} for \code{\link{flexdog}} EM algorithm.
#'
#' The key idea here is choosing the pi's so that the two modes have equal probability.
#'
#' @inheritParams flexdog
#'
#' @seealso \code{\link{flexdog}} for where this is used.
#'
#' @author David Gerard
initialize_pivec <- function(ploidy, mode, model = c("ash", "flex")) {
  model <- match.arg(model)
  if (model == "flex") {
    pivec <- rep(x = 1 / (ploidy + 1), length = ploidy + 1)
  } else if ((mode <= 0) | (mode >= ploidy) | (mode %% 1 == 0)) {
    pivec <- rep(x = 1 / (ploidy + 1), length = ploidy + 1)
  } else {
    floor_mode <- floor(mode)
    ceil_mode  <- ceiling(mode)

    d <- sum(1 / (ceil_mode:ploidy - ceil_mode + 1)) /
      (sum(1 / (floor_mode - 0:floor_mode + 1)))

    second_half <- 1 / (ploidy - ceil_mode + 1 + (floor_mode + 1) * d)

    first_half <- second_half * d

    ## assertthat::are_equal(first_half, (1 - (ploidy - ceil_mode + 1) * second_half) / (floor_mode + 1))

    pivec <- c(rep(first_half, length = floor_mode + 1),
               rep(second_half, length = ploidy - ceil_mode + 1))
  }

  return(pivec)
}

