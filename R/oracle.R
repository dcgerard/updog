## Functions for oracle calculations

#' Calculates the correlation between the true genotype and an
#' oracle esitmator.
#'
#' Calculates the correlation between the oracle MAP estimator (where we have perfect
#' knowledge about the data generation process) and the true genotype. This is a useful
#' approximation when you have a lot of individuals.
#'
#' To come up with \code{dist}, you need some additional assumptions.
#' For example, if the population is in Hardy-Weinberg equilibrium and
#' the allele frequency is \code{alpha} then you could calculate
#' \code{dist} using the R code: \code{dbinom(x = 0:ploidy, size = ploidy, prob = alpha)}.
#' Alternatively, if you know the genotypes of the individual's two parents are, say,
#' \code{ref_count1} and \code{ref_count2}, then you could use the \code{\link[updog]{get_q_array}}
#' function from the updog package: \code{get_q_array(ploidy)[ref_count1 + 1, ref_count2 + 1, ]}.
#'
#' @inheritParams oracle_mis
#'
#' @references Gerard, David, Luis Felipe Ventorim Ferrao,
#' Antonio Augusto Franco Garcia, and Matthew Stephens. 2018.
#' "Harnessing Empirical Bayes and Mendelian Segregation
#' for Genotyping Autopolyploids from Messy Sequencing Data."
#' \emph{bioRxiv}. Cold Spring Harbor Laboratory. doi:10.1101/281550.
#'
#' @export
#'
#' @author David Gerard
#'
#' @return The Pearson correlation between the true genotype and the oracle etimator.
#'
#' @examples
#' ## Hardy-Weinberg population with allele-frequency of 0.75.
#' ## Moderate bias and moderate overdispersion.
#' ## See how correlation decreases as we
#' ## increase the ploidy.
#' ploidy <- 2
#' dist <- stats::dbinom(0:ploidy, ploidy, 0.75)
#' oracle_cor(n = 100, ploidy = ploidy, seq = 0.001,
#'            bias = 0.7, od = 0.01, dist = dist)
#'
#' ploidy <- 4
#' dist <- stats::dbinom(0:ploidy, ploidy, 0.75)
#' oracle_cor(n = 100, ploidy = ploidy, seq = 0.001,
#'            bias = 0.7, od = 0.01, dist = dist)
#'
#' ploidy <- 6
#' dist <- stats::dbinom(0:ploidy, ploidy, 0.75)
#' oracle_cor(n = 100, ploidy = ploidy, seq = 0.001,
#'            bias = 0.7, od = 0.01, dist = dist)
#'
oracle_cor <- function(n,
                       ploidy,
                       seq,
                       bias,
                       od,
                       dist) {
  ## Check input --------------------------------
  assertthat::are_equal(length(n), length(ploidy), length(seq),
                        length(bias), length(od), 1)
  assertthat::assert_that(n > 0)
  assertthat::assert_that(ploidy > 0)
  assertthat::assert_that(seq >= 0, seq <= 1)
  assertthat::assert_that(bias > 0)
  assertthat::assert_that(bias >= 0, bias <= 1)
  assertthat::are_equal(length(dist), ploidy + 1)
  assertthat::assert_that(all(dist >= 0))
  assertthat::assert_that(abs(sum(dist) - 1) < .Machine$double.eps * 100)

  ## Get joint distribution --------------------
  jd <- oracle_joint(n      = n,
                     ploidy = ploidy,
                     seq    = seq,
                     bias   = bias,
                     od     = od,
                     dist   = dist)

  ## Get correlation --------------------------
  pos_geno     <- 0:ploidy
  mu_true      <- sum(dist * pos_geno)
  sigma_true   <- sqrt(sum(((pos_geno - mu_true) ^ 2) * dist))
  marge_est    <- rowSums(jd)
  mu_est       <- sum(marge_est * pos_geno)
  sigma_est    <- sqrt(sum(((pos_geno - mu_est) ^ 2) * marge_est))
  cor_true_est <- sum(jd * outer(X = pos_geno - mu_est,
                                 Y = pos_geno - mu_true,
                                 FUN = "*")) /
    (sigma_true * sigma_est)

  return(cor_true_est)
}
