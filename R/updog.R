#' \code{updog} Flexible Genotyping for Polyploids
#'
#' Implements empirical Bayes approaches to genotype
#' polyploids from next generation sequencing data while
#' accounting for allelic bias, overdispersion, and sequencing
#' error. The main function is \code{\link{flexdog}}, which allows
#' the specification of many different genotype distributions.
#' An experimental function that takes into account varying
#' levels of relatedness is implemented in \code{\link{mupdog}}.
#' Also provided are functions to simulate genotypes
#' (\code{\link{rgeno}}) and read-counts
#' (\code{\link{rflexdog}}), as well as functions to calculate
#' oracle genotyping error rates (\code{\link{oracle_mis}}) and
#' correlation with the true genotypes (\code{\link{oracle_cor}}).
#' These latter two functions are useful for read depth calculations.
#' Run \code{browseVignettes(package = "updog")} in R
#' for example usage.
#'
#' The package is named \code{updog} for "Using
#' Parental Data for Offspring Genotyping" because
#' we originally developed the
#' method for full-sib populations, but it works
#' now for more general populations.
#'
#' Our best competitor is probably the \code{fitPoly} package,
#' which you can check out at
#' \url{https://cran.r-project.org/package=fitPoly}. Though, we think
#' that \code{updog} returns better calibrated measures of uncertainty
#' when you have next-generation sequencing data.
#'
#' If you find a bug or want an enhancement, please submit an
#' issue at \url{http://github.com/dcgerard/updog/issues}.
#'
#' @references Gerard, David, Luis Felipe Ventorim Ferrao,
#' Antonio Augusto Franco Garcia, and Matthew Stephens. 2018.
#' "Harnessing Empirical Bayes and Mendelian Segregation
#' for Genotyping Autopolyploids from Messy Sequencing Data."
#' \emph{bioRxiv}. Cold Spring Harbor Laboratory. doi:10.1101/281550.
#'
#' @section \code{updog} Functions:
#' \describe{
#'   \item{\code{\link{flexdog}}}{The main function that
#'       fits an empirical Bayes approach to genotype polyploids
#'       from next generation sequencing data.}
#'   \item{\code{\link{mupdog}}}{An experimental approach to genotype
#'       autopolyploids that accounts for varying levels of
#'       relatedness between the individuals in the sample.}
#'   \item{\code{\link{rgeno}}}{simulate the genotypes of a sample
#'       from one of the models allowed in \code{\link{flexdog}}.}
#'   \item{\code{\link{rflexdog}}}{Simulate from the
#'       \code{\link{flexdog}} model.}
#'   \item{\code{\link{plot.flexdog}}}{Plotting the output of
#'       \code{\link{flexdog}}.}
#'   \item{\code{\link{plot.mupdog}}}{Plotting the output of
#'       \code{\link{mupdog}}.}
#'   \item{\code{\link{oracle_joint}}}{The joint distribution of the true genotype
#'       and an oracle estimator.}
#'   \item{\code{\link{oracle_plot}}}{Visualize the output of \code{\link{oracle_joint}}.}
#'   \item{\code{\link{oracle_mis}}}{The oracle misclassification error rate (Bayes rate).}
#'   \item{\code{\link{oracle_cor}}}{Correlation between the true genotype and the oracle estimated genotype.}
#' }
#'
#' @section \code{updog} Datasets:
#' \describe{
#'   \item{\code{\link{snpdat}}}{A small example dataset for using
#'       \code{\link{flexdog}}.}
#'   \item{\code{\link{uitdewilligen}}}{A small example dataset
#'       for using \code{\link{mupdog}}.}
#'   \item{\code{\link{mupout}}}{The output from fitting
#'       \code{\link{mupdog}} to \code{\link{uitdewilligen}}.}
#' }
#'
#' @useDynLib updog
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppArmadillo armadillo_version
#' @importFrom ggplot2 ggplot
#' @importFrom assertthat assert_that
#' @importFrom doParallel registerDoParallel
#' @importFrom ggthemes colorblind_pal
#' @importFrom foreach %dopar%
#' @importFrom stringr str_replace
#'
#' @docType package
#' @name updog
#'
#' @author David Gerard
NULL


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
