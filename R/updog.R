#' @keywords internal
#'
#' @useDynLib updog
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppArmadillo armadillo_version
#' @importFrom ggplot2 ggplot
#' @importFrom assertthat assert_that
#' @importFrom foreach %dopar%
#' @importFrom doRNG %dorng%
#' @import doFuture
#' @importFrom foreach foreach
#' @importFrom methods as
"_PACKAGE"

#' The Beta-Binomial Distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the beta-binomial distribution when parameterized
#' by the mean \code{mu} and the overdispersion parameter \code{rho}
#' rather than the typical shape parameters.
#'
#'
#' Let \eqn{\mu} and \eqn{\rho} be the mean and overdispersion parameters.
#' Let \eqn{\alpha} and \eqn{\beta} be the usual shape parameters of
#' a beta distribution. Then we have the relation
#' \deqn{\mu = \alpha/(\alpha + \beta),}
#' and
#' \deqn{\rho = 1/(1 + \alpha + \beta).}
#' This necessarily means that
#' \deqn{\alpha = \mu (1 - \rho)/\rho,}
#' and
#' \deqn{\beta = (1 - \mu) (1 - \rho)/\rho.}
#'
#' @param x,q A vector of quantiles.
#' @param p A vector of probabilities.
#' @param n The number of observations.
#' @param size A vector of sizes.
#' @param mu Either a scalar of the mean for each observation,
#'     or a vector of means of each observation, and thus
#'     the same length as \code{x} and \code{size}. This must
#'     be between 0 and 1.
#' @param rho Either a scalar of the overdispersion parameter
#'     for each observation, or a vector of overdispersion
#'     parameters of each observation, and thus the same length as
#'     \code{x} and \code{size}. This must be between 0 and 1.
#' @param log,log_p A logical vector either of length 1 or the same
#'     length as \code{x} and \code{size}. This determines whether
#'     to return the log probabilities for all observations
#'     (in the case that its length is 1) or for
#'     each observation (in the case that
#'     its length is that of \code{x} and \code{size}).
#'
#' @name betabinom
#'
#' @author David Gerard
#'
#' @return Either a random sample (\code{rbetabinom}),
#'     the density (\code{dbetabinom}), the tail
#'     probability (\code{pbetabinom}), or the quantile
#'     (\code{qbetabinom}) of the beta-binomial distribution.
#'
#' @examples
#' x <- rbetabinom(n = 10, size = 10, mu = 0.1, rho = 0.01)
#' dbetabinom(x = 1, size = 10, mu = 0.1, rho = 0.01, log = FALSE)
#' pbetabinom(q = 1, size = 10, mu = 0.1, rho = 0.01, log_p = FALSE)
#' qbetabinom(p = 0.6, size = 10, mu = 0.1, rho = 0.01)
#'
#'
NULL
