## Functions for oracle calculations

#' Calculates the correlation between the true genotype and an
#' oracle estimator.
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
#' @references
#' \itemize{
#'   \item{Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018). Genotyping Polyploids from Messy Sequencing Data. \emph{Genetics}, 210(3), 789-807. doi: \href{https://doi.org/10.1534/genetics.118.301468}{10.1534/genetics.118.301468}.}
#' }
#'
#' @export
#'
#' @author David Gerard
#'
#' @return The Pearson correlation between the true genotype and the oracle estimator.
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
  cor_true_est <- oracle_cor_from_joint(jd = jd)

  return(cor_true_est)
}


#' Calculate the correlation of the oracle estimator with the true
#' genotype from the joint distribution matrix.
#'
#' Calculates the correlation between the oracle MAP estimator (where we have perfect
#' knowledge about the data generation process) and the true genotype. This is a useful
#' approximation when you have a lot of individuals.
#'
#' @param jd A matrix of numerics. Element (i, j) is the probability
#'     of genotype i - 1 and estimated genotype j - 1. This is usually
#'     obtained from \code{\link{oracle_joint}}.
#'
#' @inherit oracle_cor return
#'
#' @seealso \code{\link{oracle_joint}} for getting \code{jd}.
#'     \code{\link{oracle_cor}} for not having to first calculate
#'     \code{jd}.
#'
#' @export
#'
#' @references
#' \itemize{
#'   \item{Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018). Genotyping Polyploids from Messy Sequencing Data. \emph{Genetics}, 210(3), 789-807. doi: \href{https://doi.org/10.1534/genetics.118.301468}{10.1534/genetics.118.301468}.}
#' }
#'
#' @author David Gerard
#'
#' @examples
#' ## Hardy-Weinberg population with allele-frequency of 0.75.
#' ## Moderate bias and moderate overdispersion.
#' ploidy <- 6
#' dist <- stats::dbinom(0:ploidy, ploidy, 0.75)
#' jd <- oracle_joint(n = 100, ploidy = ploidy, seq = 0.001,
#'                    bias = 0.7, od = 0.01, dist = dist)
#' oracle_cor_from_joint(jd = jd)
#'
#' ## Compare to oracle_cor
#' oracle_cor(n = 100, ploidy = ploidy, seq = 0.001,
#'            bias = 0.7, od = 0.01, dist = dist)
#'
#'
oracle_cor_from_joint <- function(jd) {
  assertthat::assert_that(is.matrix(jd))
  assertthat::are_equal(nrow(jd), ncol(jd))
  assertthat::assert_that(abs(sum(jd) - 1) < .Machine$double.eps * 1000)

  ploidy <- ncol(jd) - 1

  pos_geno     <- 0:ploidy
  dist         <- colSums(jd)
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


#' Get the oracle misclassification error rate directly from the
#' joint distribution of the genotype and the oracle estimator.
#'
#' @inherit oracle_mis return
#'
#' @inheritParams oracle_cor_from_joint
#'
#' @author David Gerard
#'
#' @export
#'
#' @references
#' \itemize{
#'   \item{Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018). Genotyping Polyploids from Messy Sequencing Data. \emph{Genetics}, 210(3), 789-807. doi: \href{https://doi.org/10.1534/genetics.118.301468}{10.1534/genetics.118.301468}.}
#' }
#'
#' @seealso \code{\link{oracle_joint}} for getting \code{jd}.
#'     \code{\link{oracle_mis}} for not having to first calculate
#'     \code{jd}.
#'
#' @examples
#' ## Hardy-Weinberg population with allele-frequency of 0.75.
#' ## Moderate bias and moderate overdispersion.
#' ploidy <- 6
#' dist <- stats::dbinom(0:ploidy, ploidy, 0.75)
#' jd <- oracle_joint(n = 100, ploidy = ploidy, seq = 0.001,
#'                    bias = 0.7, od = 0.01, dist = dist)
#' oracle_mis_from_joint(jd = jd)
#'
#' ## Compare to oracle_cor
#' oracle_mis(n = 100, ploidy = ploidy, seq = 0.001,
#'            bias = 0.7, od = 0.01, dist = dist)
#'
oracle_mis_from_joint <- function(jd) {
  assertthat::assert_that(is.matrix(jd))
  assertthat::are_equal(nrow(jd), ncol(jd))
  assertthat::assert_that(abs(sum(jd) - 1) < .Machine$double.eps * 1000)

  return(1 - sum(diag(jd)))
}

#' Get the oracle misclassification error rates (conditional on
#' true genotype) directly from the
#' joint distribution of the genotype and the oracle estimator.
#'
#' @inheritParams oracle_cor_from_joint
#'
#' @inherit oracle_mis_vec return
#'
#' @author David Gerard
#'
#' @export
#'
#' @references
#' \itemize{
#'   \item{Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018). Genotyping Polyploids from Messy Sequencing Data. \emph{Genetics}, 210(3), 789-807. doi: \href{https://doi.org/10.1534/genetics.118.301468}{10.1534/genetics.118.301468}.}
#' }
#'
#' @seealso \code{\link{oracle_joint}} for getting \code{jd}.
#'     \code{\link{oracle_mis_vec}} for not having to first calculate
#'     \code{jd}.
#'
#' @examples
#' ## Hardy-Weinberg population with allele-frequency of 0.75.
#' ## Moderate bias and moderate overdispersion.
#' ploidy <- 6
#' dist <- stats::dbinom(0:ploidy, ploidy, 0.75)
#' jd <- oracle_joint(n = 100, ploidy = ploidy, seq = 0.001,
#'                    bias = 0.7, od = 0.01, dist = dist)
#' oracle_mis_vec_from_joint(jd = jd)
#'
#' ## Compare to oracle_cor
#' oracle_mis_vec(n = 100, ploidy = ploidy, seq = 0.001,
#'                bias = 0.7, od = 0.01, dist = dist)
#'
oracle_mis_vec_from_joint <- function(jd) {
  assertthat::assert_that(is.matrix(jd))
  assertthat::are_equal(nrow(jd), ncol(jd))
  assertthat::assert_that(abs(sum(jd) - 1) < .Machine$double.eps * 1000)

  return(1 - diag(jd) / colSums(jd))
}



#' Construct an oracle plot from the output of \code{\link{oracle_joint}}.
#'
#' After obtaining the joint distribution of the true genotype with the estimated genotype from
#' the oracle estimator using \code{\link{oracle_joint}}, you can use \code{oracle_plot} to
#' visualize this joint distribution.
#'
#' @param jd A matrix containing the joint distribution of the true genotype and
#'     the oracle estimator. Usually, this is obtained by a call from \code{\link{oracle_joint}}.
#'
#' @author David Gerard
#'
#' @return A \code{\link[ggplot2]{ggplot}} object containing the oracle plot. The x-axis indexes
#'     the possible values of the estimated genotype. The y-axis indexes the possible values of
#'     the true genotype. The number in cell (i, j) is the probability that an individual will have
#'     true genotype i but is estimated to have genotype j. This is when using an oracle estimator.
#'     The cells are also color-coded by the size of the probability in each cell. At the top are
#'     listed the oracle misclassification error rate and the correlation of the true genotype
#'     with the estimated genotype. Both of these quantities may be derived from the joint distribution.
#'
#' @export
#'
#' @references
#' \itemize{
#'   \item{Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018). Genotyping Polyploids from Messy Sequencing Data. \emph{Genetics}, 210(3), 789-807. doi: \href{https://doi.org/10.1534/genetics.118.301468}{10.1534/genetics.118.301468}.}
#' }
#'
#' @seealso \code{\link{oracle_joint}} for obtaining \code{jd}.
#'
#' @examples
#' ploidy <- 6
#' dist <- stats::dbinom(0:ploidy, ploidy, 0.75)
#' jd <- oracle_joint(n = 100, ploidy = ploidy, seq = 0.001,
#'                    bias = 0.7, od = 0.01, dist = dist)
#' pl <- oracle_plot(jd = jd)
#' print(pl)
#'
oracle_plot <- function(jd) {
  mat_text <- format(round(jd, digits = 2), digits = 2)
  probability_text <- c(sub("0\\.", "\\.", sub("0.00", "0", mat_text)))

  dfdat <- cbind(expand.grid(x = seq_len(nrow(jd)) - 1, y = seq_len(nrow(jd)) - 1),
        probability_text,
        Probability = c(jd))

  omiss <- oracle_mis_from_joint(jd = jd)
  ocorr <- oracle_cor_from_joint(jd = jd)

  ggplot2::ggplot(data = dfdat, mapping = ggplot2::aes_string(x = "x", y = "y", fill = "Probability")) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes_string(label = "probability_text")) +
    ggplot2::xlab("Estimated Genotype") +
    ggplot2::ylab("True Genotype") +
    ggplot2::ggtitle(paste0("Misclassification Error Rate: ",
                            round(omiss, digits = 3),
                            "\nCorrelation: ",
                            round(ocorr, digits = 3))) +
    ggplot2::scale_fill_gradient2() +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = 0:(ncol(jd) - 1)) +
    ggplot2::scale_y_continuous(breaks = 0:(nrow(jd) - 1)) ->
    pl
  return(pl)
}
