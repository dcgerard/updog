################
## Code for Preferential Pairing
################

#' Probability of sending a dosage to an offspring
#'
#' This only works for tetraploids. The hope is to generalize it some day.
#'
#' This function comes from my unpublished tepp package. I might
#' remove it if I ever submit tepp to CRAN.
#'
#' @param tau The double reduciton parameter. Should be a proportion
#'     between 0 and 1 (inclusive). Should not be \code{NULL} if \code{ell} is
#'     \code{1}, \code{2}, or \code{3}.
#' @param gamma The preferential pairing parameter. Should be a proportion
#'     between 0 and 1 (inclusive). Should not be \code{NULL} if \code{ell} is
#'     \code{2}.
#' @param ell The dosage of the parent. Should be an integer between 0 and 4
#'     (inclusive).
#'
#' @return A vector of length 3. Element i is the probability that the
#'     parent sends i-1 A's to the child.
#'
#' @noRd
#'
#' @author David Gerard
prob_send_pp <- function(ell,
                         tau   = NULL,
                         gamma = NULL) {
  ## check input --------------------------------------------------------------
  stopifnot(ell >= 0, ell <= 4)
  stopifnot(length(ell) == 1)
  if (ell != 0 & ell != 4) {
    stopifnot(!is.null(tau))
    if (tau > -sqrt(.Machine$double.eps) & tau < 0) tau <- 0
    if (tau < 1 + sqrt(.Machine$double.eps) & tau > 1) tau <- 1
    stopifnot(tau >= 0, tau <= 1)
    stopifnot(length(tau) == 1)
  }
  if (ell == 2) {
    stopifnot(!is.null(gamma))
    if (gamma > -sqrt(.Machine$double.eps) & gamma < 0) gamma <- 0
    if (gamma < 1 + sqrt(.Machine$double.eps) & gamma > 1) gamma <- 1
    stopifnot(gamma >= 0, gamma <= 1)
    stopifnot(length(gamma) == 1)
  }

  probvec <- rep(NA, length.out = 3)
  if (ell == 0) {
    probvec <- c(1, 0, 0)
  } else if(ell == 1) {
    probvec[[1]] <- 0.5 + 0.25 * tau
    probvec[[2]] <- 0.5 - 0.5 * tau
    probvec[[3]] <- 0.25 * tau
  } else if(ell == 2) {
    probvec[[1]] <- 0.5 * tau + 0.25 * (1 - tau) * (1 - gamma)
    probvec[[2]] <- 0.5 * (1 - tau) * (1 + gamma)
    probvec[[3]] <- probvec[[1]]
  } else if(ell == 3) {
    probvec[[1]] <- 0.25 * tau
    probvec[[2]] <- 0.5 - 0.5 * tau
    probvec[[3]] <- 0.5 + 0.25 * tau
  } else if(ell == 4) {
    probvec <- c(0, 0, 1)
  }
  return(probvec)
}


#' Probability distribution of offspring dosages given parental parameters
#'
#' Only works for tetraploids.
#'
#' This function comes from my unpublished tepp package. I might
#' remove it if I ever submit tepp to CRAN.
#'
#' @param ell1 Dosage of parent 1. Should be an integer between 0 and 4
#'     (inclusive).
#' @param ell2 Dosage of parent 2. Should be an integer between 0 and 4
#'     (inclusive).
#' @param tau1 Double reduction parameter of parent 1. Should be a proportion
#'     between 0 and 1 (inclusive). Should not be \code{NULL} if \code{ell1} is
#'     \code{1}, \code{2}, or \code{3}.
#' @param tau2 Double reduction parameter of parent 2. Should be a proportion
#'     between 0 and 1 (inclusive). Should not be \code{NULL} if \code{ell2} is
#'     \code{1}, \code{2}, or \code{3}.
#' @param gamma1 Preferential pairing parameter of parent 1. Should be a proportion
#'     between 0 and 1 (inclusive). Should not be \code{NULL} if \code{ell1} is
#'     \code{2}.
#' @param gamma2 Preferential pairing parameter of parent 2. Should be a proportion
#'     between 0 and 1 (inclusive). Should not be \code{NULL} if \code{ell2} is
#'     \code{2}.
#'
#' @return A vector of length 5. Element i is the probability that a child
#'     will have dosage i-1.
#'
#' @noRd
#'
#' @author David Gerard
prob_dosage_pp <- function(ell1,
                           ell2,
                           tau1   = NULL,
                           tau2   = NULL,
                           gamma1 = NULL,
                           gamma2 = NULL) {

  ## parameters are checked in prob_send()
  p1send <- prob_send_pp(ell = ell1, tau = tau1, gamma = gamma1)
  p2send <- prob_send_pp(ell = ell2, tau = tau2, gamma = gamma2)

  ## convolve
  dosedist <- stats::convolve(x = p1send, y = rev(p2send), type = "open")

  return(dosedist)
}

#' Mixture of discrete uniform and \code{\link{prob_send_pp}()}
#'
#'
#' @inheritParams prob_dosage_pp
#' @param fs1alpha The mixing proportion. This is the probability
#'     of a discrete uniform.
#'
#' @author David Gerard
#'
#' @noRd
prob_dosage_pp_unif <- function(ell1,
                                ell2,
                                tau1     = NULL,
                                tau2     = NULL,
                                gamma1   = NULL,
                                gamma2   = NULL,
                                fs1alpha = 0.001) {
  stopifnot(fs1alpha >= 0 & fs1alpha <= 1)

  pvec <- prob_dosage_pp(ell1   = ell1,
                         ell2   = ell2,
                         tau1   = tau1,
                         tau2   = tau2,
                         gamma1 = gamma1,
                         gamma2 = gamma2)

  retvec <- pvec * (1 - fs1alpha) + fs1alpha / length(pvec)

  return(retvec)
}

#' Obtain parameters from vectorization of parameters.
#'
#' Only works in tetraploids
#'
#' @param par A vector of parameters in preference order of
#'     \code{c(tau1, gamma1, tau2, gamma2)}. But the actual elements
#'     depend on the values of \code{ell}.
#' @param ell1 Parental 1 genotype.
#' @param ell2 Parental 2 genotype.
#'
#' @author David Gerard
#'
#' @noRd
undo_par <- function(par, ell1, ell2) {
  retvec <- rep(NA_real_, length = 4)
  names(retvec) <- c("tau1", "tau2", "gamma1", "gamma2")
  if(ell1 %in% c(0, 4) & ell2 %in% c(0, 4)) {
    # do nothing
  } else if (ell1 %in% c(1, 3) & ell2 %in% c(0, 4)) {
    stopifnot(length(par) == 1)
    retvec[["tau1"]] <- par[[1]]
  } else if (ell1 %in% c(0, 4) & ell2 %in% c(1, 3)) {
    stopifnot(length(par) == 1)
    retvec[["tau2"]] <- par[[1]]
  } else if (ell1 %in% c(1, 3) & ell2 %in% c(1, 3)) {
    stopifnot(length(par) == 2)
    retvec[["tau1"]] <- par[[1]]
    retvec[["tau2"]] <- par[[2]]
  } else if (ell1 == 2 & ell2 %in% c(0, 4)) {
    stopifnot(length(par) == 2)
    retvec[["tau1"]] <- par[[1]]
    retvec[["gamma1"]] <- par[[2]]
  } else if (ell1 %in% c(0, 4) & ell2 == 2) {
    stopifnot(length(par) == 2)
    retvec[["tau2"]] <- par[[1]]
    retvec[["gamma2"]] <- par[[2]]
  } else if (ell1 == 2 & ell2 %in% c(1, 3)) {
    stopifnot(length(par) == 3)
    retvec[["tau1"]] <- par[[1]]
    retvec[["gamma1"]] <- par[[2]]
    retvec[["tau2"]] <- par[[3]]
  } else if (ell1 %in% c(1, 3) & ell2 == 2) {
    stopifnot(length(par) == 3)
    retvec[["tau1"]] <- par[[1]]
    retvec[["tau2"]] <- par[[2]]
    retvec[["gamma2"]] <- par[[3]]
  } else if (ell1 == 2 & ell2 == 2) {
    stopifnot(length(par) == 4)
    retvec[["tau1"]] <- par[[1]]
    retvec[["gamma1"]] <- par[[2]]
    retvec[["tau2"]] <- par[[3]]
    retvec[["gamma2"]] <- par[[4]]
  }
  return(retvec)
}


#' Sets vectorization formation given parameters.
#'
#' @param ell1 Parental genotype 1
#' @param ell2 Parental genotype 2
#' @param tau1 The double reduction parameter for parent 1
#' @param tau2 The double reduction parameter for parent 2
#' @param gamma1 The preferential pairing parameter for parent 1.
#' @param gamma2 The preferential pairing parameter for parent 2.
#' @param tau_lower The lower bound on tau1 and tau2
#' @param tau_upper The upper bound on tau1 and tau2
#' @param gamma1_lower The lower bound on gamma1 and gamma2
#' @param gamma2_upper The upper bound on gamma1 and gamma2
#'
#' @author David Gerard
#'
#' @return A matrix. The first row is the vector par. The second row
#'     is the vector of lower bounds. The third row is the vector of
#'     upper bounds
#'
#' @noRd
set_par <- function(ell1,
                    ell2,
                    tau1 = NULL,
                    tau2 = NULL,
                    gamma1 = NULL,
                    gamma2 = NULL,
                    tau_lower = 0,
                    tau_upper = 1/6,
                    gamma_lower = 0,
                    gamma_upper = 1) {
  if(ell1 %in% c(0, 4) & ell2 %in% c(0, 4)) {
    par <- NULL
    lower <- NULL
    upper <- NULL
  } else if (ell1 %in% c(1, 3) & ell2 %in% c(0, 4)) {
    par <- c(tau1)
    names(par) <- "tau1"
    lower <- c(tau_lower)
    upper <- c(tau_upper)
  } else if (ell1 %in% c(0, 4) & ell2 %in% c(1, 3)) {
    par <- c(tau2)
    names(par) <- "tau2"
    lower <- c(tau_lower)
    upper <- c(tau_upper)
  } else if (ell1 %in% c(1, 3) & ell2 %in% c(1, 3)) {
    par <- c(tau1, tau2)
    names(par) <- c("tau1", "tau2")
    lower <- c(tau_lower, tau_lower)
    upper <- c(tau_upper, tau_upper)
  } else if (ell1 == 2 & ell2 %in% c(0, 4)) {
    par <- c(tau1, gamma1)
    names(par) <- c("tau1", "gamma1")
    lower <- c(tau_lower, gamma_lower)
    upper <- c(tau_upper, gamma_upper)
  } else if (ell1 %in% c(0, 4) & ell2 == 2) {
    par <- c(tau2, gamma2)
    names(par) <- c("tau2", "gamma2")
    lower <- c(tau_lower, gamma_lower)
    upper <- c(tau_upper, gamma_upper)
  } else if (ell1 == 2 & ell2 %in% c(1, 3)) {
    par <- c(tau1, gamma1, tau2)
    names(par) <- c("tau1", "gamma1", "tau2")
    lower <- c(tau_lower, gamma_lower, tau_lower)
    upper <- c(tau_upper, gamma_upper, tau_upper)
  } else if (ell1 %in% c(1, 3) & ell2 == 2) {
    par <- c(tau1, tau2, gamma2)
    names(par) <- c("tau1", "tau2", "gamma2")
    lower <- c(tau_lower, tau_lower, gamma_lower)
    upper <- c(tau_upper, tau_upper, gamma_upper)
  } else if (ell1 == 2 & ell2 == 2) {
    par <- c(tau1, gamma1, tau2, gamma2)
    names(par) <- c("tau1", "gamma1", "tau2", "gamma2")
    lower <- c(tau_lower, gamma_lower, tau_lower, gamma_lower)
    upper <- c(tau_upper, gamma_upper, tau_upper, gamma_upper)
  }
  return(rbind(par, lower, upper))
}

#' Objective function that is optimized in \code{\link{update_f1_s1_pp}()}
#'
#' This is \deqn{\sum_{k=0}^Kw_f(k|tau1, gamma1, tau2, gamma2, ell1, ell2)}
#' for given weights \eqn{w_k}.
#'
#' @param par The parameter vector. See \code{\link{undo_par}()} for a description.
#' @param weight_vec The \eqn{w_k}'s.
#' @param ell1 The parent 1 dosage.
#' @param ell2 The parent 2 dosage
#' @param fs1alpha The uniform mixing weight to avoid degeneracy.
#' @param p1pen This should be log-BB(x|ell1, seq, bias, od) for parent 1 for each ell1.
#'     So it's a vector.
#' @param p2pen This should be log-BB(x|ell1, seq, bias, od) for parent 2 for each ell2.
#'     So it's a vector.
#' @param pop Either F1 or S1. We copy the parameters if s1.
#' @author David Gerard
#'
#' @noRd
f1_s1_pp_obj <- function(par,
                         weightvec,
                         ell1,
                         ell2,
                         fs1alpha = 0.001,
                         p1pen = rep(0, length(weightvec)),
                         p2pen = rep(0, length(weightvec)),
                         pop = c("f1", "s1")) {
  if (is.null(p1pen)) {
    p1pen <- rep(0, length(weightvec))
  }
  if (is.null(p2pen)) {
    p2pen <- rep(0, length(weightvec))
  }
  stopifnot(length(weightvec) == length(p1pen),
            length(weightvec) == length(p2pen))

  if (pop == "s1") { ## expect par to be half as long if s1
    par <- c(par, par)
  }

  paramvec <- undo_par(par = par, ell1 = ell1, ell2 = ell2)
  priorvec <- prob_dosage_pp_unif(ell1 = ell1,
                                  ell2 = ell2,
                                  tau1 = paramvec[["tau1"]],
                                  tau2 = paramvec[["tau2"]],
                                  gamma1 = paramvec[["gamma1"]],
                                  gamma2 = paramvec[["gamma2"]],
                                  fs1alpha = fs1alpha)
  sum(log(priorvec) * weightvec) + p1pen[[ell1 + 1]] + p2pen[[ell2 + 1]]
}

#' Maximize the objective function in \code{\link{f1_s1_pp_obj}()}
#'
#' The objective function is
#' \deqn{\sum_{k=0}^Kw_f(k|tau1, gamma1, tau2, gamma2, ell1, ell2)}
#' for given weights \eqn{w_k}.
#'
#' Does a bunch of optim() calls. Fewer if pop = "s1".
#'
#' @param weight_vec The \eqn{w_k}'s
#' @param pop Is this an F1 or S1 population?
#' @param tau1_init The initialization for tau1.
#' @param tau2_init The initialization for tau2
#' @param gamma1_init The initialization for gamma1
#' @param gamma2_init The initialization for gamma2
#' @param fs1alpha The uniform mixing weight to avoid degeneracy
#' @param p1pen This should be log-BB(x|ell1, seq, bias, od) for parent 1.
#' @param p2pen This should be log-BB(x|ell1, seq, bias, od) for parent 2.
#'
#' @author David Gerard
#'
#' @noRd
update_f1_s1_pp <- function(weightvec,
                            pop = c("f1", "s1"),
                            tau1_init = 0.05,
                            tau2_init = 0.05,
                            gamma1_init = 1/3,
                            gamma2_init = 1/3,
                            fs1alpha = 0.001,
                            p1pen = rep(0, length(weightvec)),
                            p2pen = rep(0, length(weightvec))) {
  pop <- match.arg(pop)
  ploidy <- 4

  obest <- list()
  obest$value <- -Inf
  for (ell1 in 0:ploidy){
    if (pop == "s1") {
      ell2_vec <- ell1
    } else if (pop == "f1") {
      ell2_vec <- 0:ploidy
    } else {
      stop("how did you get here?")
    }

    for (ell2 in ell2_vec) {
      if (ell1 %in% c(0, 4) & ell2 %in% c(0, 4)) {
        oout <- list()
        oout$value <- f1_s1_pp_obj(par = NULL,
                                   weightvec = weightvec,
                                   ell1 = ell1,
                                   ell2 = ell2,
                                   fs1alpha = fs1alpha,
                                   p1pen = p1pen,
                                   p2pen = p2pen,
                                   pop = pop)
      } else {
        parmat <- set_par(ell1        = ell1,
                          ell2        = ell2,
                          tau1        = tau1_init,
                          tau2        = tau2_init,
                          gamma1      = gamma1_init,
                          gamma2      = gamma2_init,
                          tau_lower   = 0,
                          tau_upper   = 1/6,
                          gamma_lower = 0,
                          gamma_upper = 1)

        if (pop == "s1") {
          parmat <- parmat[, seq_len(ncol(parmat) / 2), drop = FALSE]
        }

        oout <- stats::optim(par       = parmat[1, ],
                             fn        = f1_s1_pp_obj,
                             method    = "L-BFGS-B",
                             lower     = parmat[2, ],
                             upper     = parmat[3, ],
                             control   = list(fnscale = -1, maxit = 10),
                             weightvec = weightvec,
                             ell1      = ell1,
                             ell2      = ell2,
                             fs1alpha  = fs1alpha,
                             p1pen     = p1pen,
                             p2pen     = p2pen,
                             pop       = pop)
      }
      oout$ell1 <- ell1
      oout$ell2 <- ell2
      if (oout$value > obest$value) {
        obest <- oout
      }
    }
  }

  if (pop == "f1") {
    retvec <- c(undo_par(par = obest$par, ell1 = obest$ell1, ell2 = obest$ell2),
                ell1 = obest$ell1, ell2 = obest$ell2)
  } else if (pop == "s1") {
    stopifnot(obest$ell1 == obest$ell2)
    retvec <- undo_par(par = c(obest$par, obest$par), ell1 = obest$ell1, obest$ell2)
    retvec <- c(retvec[names(retvec) %in% c("tau1", "gamma1")], ell1 = obest$ell1)
  } else {
    stop("how did you get here?")
  }

  return(retvec)
}
