

#' KL-divergence between a provided distribution and a Gaussian-Binomial distribution.
#'
#' @param dist A vector of probabilities. The \eqn{i}th element is the probability of having
#'     \eqn{i+1} copies of the reference allele.
#' @param ploidy The ploidy of the individual.
#' @param sigma2 The underlying variance of the Gaussian-Binomial distribution.
#' @param alpha The probability of success.
#' 
#' @return A double. The KL-divergence between \code{dist} and the
#'     Guassian-Binomial. 
#'
#' @author David Gerard
#'
#' @seealso \code{\link{kl_opt}} for optimizing \code{kl_prob} over both \code{sigma2} and \code{alpha}.
#'
#'     \code{\link{kl_wrapp}} for a wrapper of \code{kl_prob}.
#'
kl_prob <- function(dist, ploidy, sigma2, alpha) {
  assertthat::are_equal(length(dist), ploidy + 1)
  assertthat::are_equal(length(alpha), 1)
  assertthat::are_equal(length(sigma2), 1)
  assertthat::are_equal(length(ploidy), 1)
  assertthat::assert_that(sigma2 > 0)
  assertthat::assert_that(alpha > 0, alpha < 1)
  assertthat::assert_that(ploidy > 0)
  assertthat::are_equal(ploidy %% 1, 0)
  assertthat::are_equal(sum(dist), 1)
  assertthat::assert_that(all(dist >= 0))
  assertthat::assert_that(all(dist <= 1))

  mu <- matrix(data = 0, nrow = 1, ncol = 1)
  rho <- 0
  sigma2 <- matrix(data = sigma2, nrow = 1, ncol = 1)
  prob_vec <- c(compute_all_post_prob(ploidy = ploidy, mu = mu, sigma2 = sigma2, alpha = alpha, rho = rho))
  if (any(prob_vec == 0)) {
    prob_vec[prob_vec == 0] <- .Machine$double.eps
  }
  return(-1 * sum(dist * log(prob_vec)))
}

#' A wrapper of \code{\link{kl_prob}}.
#'
#' @inheritParams kl_prob
#' @inherit kl_prob return
#' @param par A vector of length two. The first element is \code{sigma2} in \code{\link{kl_prob}},
#'     the second element is \code{alpha} in \code{\link{kl_prob}}.
#'
#' @author David Gerard
#' 
#'
#' @seealso \code{\link{kl_prob}} and \code{\link{kl_opt}} for optimizing \code{kl_wrapp} over \code{par}.
#'
kl_wrapp <- function(par, dist, ploidy) {
  assertthat::are_equal(length(par), 2)
  kl_prob(dist = dist, ploidy = ploidy, sigma2 = par[1], alpha = par[2])
}

#' Finds the Gaussian-Binomial distribution closest (in terms of Kullback-Leibler divergence)
#' to a provided distribution.
#'
#' The Gaussian-Binomial distribution is defined here by the hierarchical model
#' \deqn{z ~ N(0, \sigma^2),}
#' \deqn{y = F^{-1}(\Phi(z)|K, \alpha),}
#' where \eqn{F^{-1}(.|K, \alpha)} is the quantile function of a Binomial distribution with size
#' \eqn{K} and success probability \eqn{\alpha}. \eqn{\Phi(.)} is the standard normal cumulative
#' distribution function. In our context, \eqn{K} is the ploidy of the species. This function tries to
#' find the closest Gaussian-Binomial distribution to a provided distribution, where "closest" means
#' smallest Kulback-Leibler divergence. We do this by optimizing over both \eqn{\alpha} and \eqn{\sigma^2}.
#'
#' @return A list containing the fitted distribution along with the 
#'     fitted parameters.
#'
#' @inheritParams kl_prob
#'
#' @seealso \code{\link{kl_prob}} and \code{\link{kl_wrapp}} for the objective functions.
#'
#' @author David Gerard
#'
#'
kl_opt <- function(dist, ploidy) {
  assertthat::are_equal(length(dist), ploidy + 1)
  assertthat::are_equal(sum(dist), 1)
  alpha_start <- sum((0:ploidy) * dist) / ploidy
  sigma2_start <- 1

  tol <- .Machine$double.eps * 10

  oout <- stats::optim(par = c(sigma2_start, alpha_start), fn = kl_wrapp,
               lower = c(tol, tol), upper = c(Inf, 1 - tol), method = "L-BFGS-B",
               dist = dist, ploidy = ploidy)

  kl_wrapp(par = c(tol, tol), dist = dist, ploidy = ploidy)

  fitted_dist <- c(compute_all_post_prob(ploidy = ploidy, mu = matrix(data = 0, nrow = 1, ncol = 1),
                                         sigma2 = matrix(data = oout$par[1], nrow = 1, ncol = 1),
                                         alpha = oout$par[2], rho = 0))
  return_list <- list()
  return_list$sigma2 <- oout$par[1]
  return_list$alpha  <- oout$par[2]
  return_list$dist   <- fitted_dist
  return(return_list)
}
