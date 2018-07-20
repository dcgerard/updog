## Functions related to double reduction.

#' Double reduction version of \code{\link{get_bivalent_probs}}.
#'
#' @inherit get_bivalent_probs
#'
#' @return A list. The same elements as in \code{\link{get_bivalent_probs}}, augmented to include
#'     more scenarios,
#'     but with the additional element \code{penvec}. This is a logical vector that
#'     is \code{TRUE} if the corresponding rows of \code{probmat} and \code{pmat} and
#'     elements of \code{lvec} would be
#'     included in the non-double-reduction-model and is \code{FALSE} otherwise.
#'
#' @seealso \code{\link{get_bivalent_probs}}.
#'
#' @author David Gerard
get_bivalent_probs_dr <- function(ploidy) {
  blist <- get_bivalent_probs(ploidy = ploidy)

  contrib_mat <- matrix(NA, ncol = ploidy + 1, nrow = nrow(blist$pmat))
  dimnames(contrib_mat) <- list(configuration = rownames(blist$probmat),
                                ell = 0:ploidy)
  probmat <- matrix(NA, nrow = 0, ncol = ploidy / 2 + 1)
  pmat    <- matrix(NA, nrow = 0, ncol = 3)
  lvec    <- vector(mode = "numeric", length = 0)
  penvec  <- vector(mode = "logical", length = 0)
  for (ell in 0:ploidy) {
    if (ell <= ploidy / 2) {
      contrib_vec <- blist$lvec <= 2 * ell
    } else {
      contrib_vec <- ploidy - blist$lvec <= 2 * (ploidy - ell)
    }
    probmat <- rbind(probmat, blist$probmat[contrib_vec, , drop = FALSE])
    pmat <- rbind(pmat, blist$pmat[contrib_vec, , drop = FALSE])
    lvec <- c(lvec, rep(ell, length = sum(contrib_vec)))
    penvec <- c(penvec, rowSums(sweep(x = blist$pmat[contrib_vec, , drop = FALSE],
                                      MARGIN = 2,
                                      STATS = c(0, 1, 2),
                                      FUN = `*`)) == ell)
  }
  dblist <- list(probmat = probmat, pmat = pmat, lvec = lvec, penvec = penvec)
  return(dblist)
}

#' Penalty used in \code{\link{update_dr}}.
#'
#' A dirichlet prior on \code{pairweights}. Returns log density.
#'
#' @param pairweights The mixing proportions to penalize.
#' @param mixing_pen The corresponding penalties.
#'
#' @seealso \code{\link{update_dr}}
#'
#' @author David Gerard
dr_pen <- function(pairweights, mixing_pen) {
  sum((log(pairweights) * mixing_pen)[mixing_pen > 10^-6])
}


#' Same as \code{\link{update_pp_f1}} but I exclusively use the EM (instead of also Brent's method),
#' and I allow for priors on the mixing proportions.
#'
#' @inherit update_pp_f1
#' @param model The model to assume.
#'
#' @seealso \code{\link{update_pp_f1}}
update_dr <- function(weight_vec,
                      model = c("f1pp", "f1ppdr"),
                      control) {
  model <- match.arg(model)
  model <- match.arg(model)
  assertthat::assert_that(!is.null(control$blist))
  assertthat::assert_that(!is.null(control$fs1_alpha))
  assertthat::assert_that(!is.null(control$p1_pair_weights))
  assertthat::assert_that(!is.null(control$pivec))
  if (model == "f1pp" | model == "f1ppdr") {
    stopifnot(!is.null(control$p2_pair_weights))
  } else {
    control$p2_pair_weights <- control$p1_pair_weights
  }

  ploidy <- length(weight_vec) - 1

  ## number of mixing components for each genotype
  num_comp <- table(control$blist$lvec)
  return_list <- list()
  return_list$p1_pair_weights <- control$p1_pair_weights
  return_list$p2_pair_weights <- control$p2_pair_weights
  return_list$obj             <- -Inf
  return_list$p1geno          <- NA
  return_list$p2geno          <- NA
  return_list$pivec           <- rep(NA, length = ploidy + 1)

  for (p1geno in 0:ploidy) {
    if ((model == "f1pp" | model == "f1ppdr") & (!is.null(control$p1_lbb) | !is.null(control$p2_lbb))) {
      p2geno_vec <- 0:ploidy
    } else if (model == "f1pp" | model == "f1ppdr") {
      p2geno_vec <- p1geno:ploidy
    } else {
      p2geno_vec <- p1geno
    }
    ## for pivec, I don't mix with uniform until very end ----
    for(p2geno in p2geno_vec) {
      temp_list <- list()
      if ((num_comp[p1geno + 1] == 1) & (num_comp[p2geno + 1] == 1)) {
        ## Just return likelihood
        p1_seg_prob <- control$blist$probmat[control$blist$lvec == p1geno, , drop = TRUE]
        p2_seg_prob <- control$blist$probmat[control$blist$lvec == p2geno, , drop = TRUE]
        temp_list$pivec <- c(convolve(p1_seg_prob, p2_seg_prob))
        temp_list$obj <- f1_obj(alpha = control$fs1_alpha,
                                pvec = temp_list$pivec,
                                weight_vec = weight_vec)
      } else if (num_comp[p1geno + 1] == 1) {
        ## EM on 2
        p1_seg_prob <- control$blist$probmat[control$blist$lvec == p1geno, , drop = TRUE]
        p2_seg_prob_mat <- control$blist$probmat[control$blist$lvec == p2geno, , drop = FALSE]

        lmat <- get_conv_inner_weights(psegprob = p1_seg_prob,
                                       psegmat  = p2_seg_prob_mat)

        initial_comp_weights <- control$p2_pair_weights[[p2geno + 1]]

        return_list$p2_pair_weights[[p2geno + 1]] <-
          c(uni_em_const(weight_vec = weight_vec,
                         lmat       = lmat,
                         pi_init    = initial_comp_weights,
                         alpha      = control$fs1_alpha,
                         lambda     = control$mixing_pen[control$blist$lvec == p2geno],
                         itermax    = 100,
                         obj_tol    = 10^-5))

        temp_list$pivec <- pivec_from_segmats(p1segmat = matrix(p1_seg_prob, nrow = 1),
                                              p2segmat = p2_seg_prob_mat,
                                              p1weight = 1,
                                              p2weight = return_list$p2_pair_weights[[p2geno + 1]])
        temp_list$obj <- uni_obj_const(pivec = return_list$p2_pair_weights[[p2geno + 1]],
                                       alpha = control$fs1_alpha,
                                       weight_vec = weight_vec,
                                       lmat = lmat,
                                       lambda = control$mixing_pen[control$blist$lvec == p2geno])

      } else if (num_comp[p2geno + 1] == 1) {
        ## EM on 1
        p2_seg_prob <- control$blist$probmat[control$blist$lvec == p2geno, , drop = TRUE]
        p1_seg_prob_mat <- control$blist$probmat[control$blist$lvec == p1geno, , drop = FALSE]

        lmat <- get_conv_inner_weights(psegprob = p2_seg_prob,
                                       psegmat  = p1_seg_prob_mat)

        initial_comp_weights <- control$p1_pair_weights[[p1geno + 1]]

        return_list$p1_pair_weights[[p1geno + 1]] <-
          c(uni_em_const(weight_vec = weight_vec,
                         lmat       = lmat,
                         pi_init    = initial_comp_weights,
                         alpha      = control$fs1_alpha,
                         lambda     = control$mixing_pen[control$blist$lvec == p1geno],
                         itermax    = 100,
                         obj_tol    = 10^-5))

        temp_list$pivec <- pivec_from_segmats(p2segmat = matrix(p2_seg_prob, nrow = 1),
                                              p1segmat = p1_seg_prob_mat,
                                              p2weight = 1,
                                              p1weight = return_list$p1_pair_weights[[p1geno + 1]])
        temp_list$obj <- uni_obj_const(pivec = return_list$p1_pair_weights[[p1geno + 1]],
                                       alpha = control$fs1_alpha,
                                       weight_vec = weight_vec,
                                       lmat = lmat,
                                       lambda = control$mixing_pen[control$blist$lvec == p1geno])
      } else {
        ## EM on both
        p1_seg_prob_mat   <- control$blist$probmat[control$blist$lvec == p1geno, , drop = FALSE]
        p2_seg_prob_mat   <- control$blist$probmat[control$blist$lvec == p2geno, , drop = FALSE]
        p1_weight_vec     <- control$p1_pair_weights[[p1geno + 1]]
        p2_weight_vec     <- control$p2_pair_weights[[p2geno + 1]]

        brent_obj     <- -Inf
        brent_tol     <- 10^-5
        brent_err     <- brent_tol + 1
        brent_iter    <- 1
        brent_itermax <- 100
        while ((brent_err > brent_tol) & (brent_iter < brent_itermax)) {
          brent_obj_old <- brent_obj

          # Update p1_weight_vec
          p2_seg_prob <- colSums(p2_weight_vec * p2_seg_prob_mat)

          lmat <- get_conv_inner_weights(psegprob = p2_seg_prob,
                                         psegmat  = p1_seg_prob_mat)

          p1_weight_vec <- c(uni_em_const(weight_vec = weight_vec,
                                          lmat       = lmat,
                                          pi_init    = p1_weight_vec,
                                          alpha      = control$fs1_alpha,
                                          lambda     = control$mixing_pen[control$blist$lvec == p1geno],
                                          itermax    = 100,
                                          obj_tol    = 10^-5))
          # Update p2_weight_vec
          p1_seg_prob <- colSums(p1_weight_vec * p1_seg_prob_mat)

          lmat <- get_conv_inner_weights(psegprob = p1_seg_prob,
                                         psegmat  = p2_seg_prob_mat)

          p2_weight_vec <- c(uni_em_const(weight_vec = weight_vec,
                                          lmat       = lmat,
                                          pi_init    = p2_weight_vec,
                                          alpha      = control$fs1_alpha,
                                          lambda     = control$mixing_pen[control$blist$lvec == p2geno],
                                          itermax    = 100,
                                          obj_tol    = 10^-5))


          brent_obj <- uni_obj_const(pivec = p2_weight_vec,
                                     alpha = control$fs1_alpha,
                                     weight_vec = weight_vec,
                                     lmat = lmat,
                                     lambda = control$mixing_pen[control$blist$lvec == p2geno])
          brent_obj <- brent_obj + dr_pen(pairweights = p1_weight_vec, mixing_pen = control$mixing_pen[control$blist$lvec == p1geno])

          ## Check same as brent_obj
          # uni_obj_const(pivec = p2_weight_vec,
          #               alpha = control$fs1_alpha,
          #               weight_vec = weight_vec,
          #               lmat = lmat,
          #               lambda = 0) +
          #   dr_pen(pairweights = p1_weight_vec, mixing_pen = control$mixing_pen[control$blist$lvec == p1geno]) +
          #   dr_pen(pairweights = p2_weight_vec, mixing_pen = control$mixing_pen[control$blist$lvec == p2geno])


          brent_err <- abs(brent_obj - brent_obj_old)
          if (brent_obj < brent_obj_old - 10^-5) {
            stop("update_dr: objective not increasing")
          }
          brent_iter <- brent_iter + 1
        }

        temp_list$obj <- brent_obj

        return_list$p1_pair_weights[[p1geno + 1]] <- p1_weight_vec
        return_list$p2_pair_weights[[p2geno + 1]] <- p2_weight_vec

        temp_list$pivec <- pivec_from_segmats(p1segmat = p1_seg_prob_mat,
                                              p2segmat = p2_seg_prob_mat,
                                              p1weight = return_list$p1_pair_weights[[p1geno + 1]],
                                              p2weight = return_list$p2_pair_weights[[p2geno + 1]])
      }
      ## Add parental contributions -------------------
      if (!is.null(control$p1_lbb)) {
        temp_list$obj <- temp_list$obj + control$p1_lbb[p1geno + 1]
      }
      if (!is.null(control$p2_lbb)) {
        temp_list$obj <- temp_list$obj + control$p2_lbb[p2geno + 1]
      }

      ## Check to see what has the highest likelihood ----
      if (temp_list$obj > return_list$obj) {
        return_list$pivec  <- temp_list$pivec
        return_list$obj    <- temp_list$obj
        return_list$p1geno <- p1geno
        return_list$p2geno <- p2geno
      }
    }
  }
  return_list$pivec <- (1 - control$fs1_alpha) * return_list$pivec +
    control$fs1_alpha / (ploidy + 1)
  return(return_list)
}
