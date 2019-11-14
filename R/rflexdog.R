## Function to simulate from the flexdog model

#' Simulate individual genotypes from one of the supported \code{\link{flexdog}} models.
#'
#' This will simulate genotypes of a sample of individuals drawn from one of the populations supported by
#' \code{\link{flexdog}}. See the details of \code{\link{flexdog}} for the models allowed. These genotype
#' distributions are described in detail in Gerard and Ferrão (2019).
#'
#' The allowable variable values of \code{allele_freq}, \code{od}, \code{p1geno}, \code{p2geno}, \code{pivec},
#' and \code{mode} varies based on the value of \code{model}. If \code{model = "ash"}, then only
#' \code{mode} and \code{pivec} can be non-\code{NULL}. If \code{model = "flex"} then only
#' \code{pivec} can be non-\code{NULL}. If \code{model = "hw"}, then only \code{allele_freq} can
#' be non-\code{NULL}. If \code{model = "f1"} then only \code{p1geno} and \code{p2geno} can be non-\code{NULL}.
#' If \code{model = "s1"} then only \code{p1geno} can be non-\code{NULL}. If \code{model = "uniform"}, then
#' none of the above variables can be non-\code{NULL}. If \code{model = "bb"}, then only \code{allele_freq},
#' and \code{od} can be non-\code{NULL}. If \code{model == "norm"}, then only \code{mu} and \code{sigma}
#' can be non-\code{NULL}.
#'
#' @param n The number of observations.
#' @param ploidy The ploidy of the species.
#' @param model What form should the prior take? See Details in \code{\link{flexdog}}.
#' @param allele_freq If \code{model = "hw"}, then this is the allele frequency of the population.
#'     For any other model, this should be \code{NULL}.
#' @param od If \code{model = "bb"}, then this is the overdispersion parameter of the beta-binomial
#'     distribution. See \code{\link{betabinom}} for details. For any other model, this should be
#'     \code{NULL}.
#' @param p1geno Either the first parent's genotype if \code{model = "f1"}
#'     (or \code{model = "f1pp"}), or the only parent's
#'     genotype if \code{model = "s1"} (or \code{model = "s1pp"}).
#'     For any other model, this should be \code{NULL}.
#' @param p2geno The second parent's genotype if \code{model = "f1"}
#'     (or \code{model = "f1pp"}).
#'     For any other model, this should be \code{NULL}.
#' @param pivec A vector of probabilities. If \code{model = "ash"}, then this represents
#'     the mixing proportions of the discrete uniforms. If
#'     \code{model = "flex"}, then element \code{i} is the probability of genotype \code{i - 1}.
#'     For any other model, this should be \code{NULL}.
#' @param mode If \code{model = "ash"}, this is the center of the distribution. This should be a non-integer value
#'     (so the mode is either the floor or the ceiling of \code{mode}).
#'     For any other model, this should be \code{NULL}.
#' @param mu If \code{model = "norm"}, this is the mean of the normal.
#'     For any other model, this should be \code{NULL}.
#' @param sigma If \code{model = "norm"}, this is the standard deviation of the normal.
#'     For any other model, this should be \code{NULL}.
#' @param p1_pair_weights The mixing weights for the bivalent
#'     pairs output in \code{\link{get_bivalent_probs}}
#'     at the \code{lvec} level of \code{p1geno}.
#'     If \code{model = "f1pp"} then this is for the first parent.
#'     If \code{model = "s1pp"}, then this is for the only parent.
#'     This should be \code{NULL} for all other values of \code{model}.
#' @param p2_pair_weights If \code{model = "s1pp"},
#'     these are the mixing weights for the bivalent
#'     pairs output in \code{\link{get_bivalent_probs}} at the
#'     \code{lvec} level of \code{p2geno} for the second
#'     parent.
#'     This should be \code{NULL} for all other values of \code{model}.
#'
#'
#' @return A vector of length \code{n} with the genotypes of the sampled individuals.
#'
#' @references
#' \itemize{
#'   \item{Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018). Genotyping Polyploids from Messy Sequencing Data. \emph{Genetics}, 210(3), 789-807. doi: \href{https://doi.org/10.1534/genetics.118.301468}{10.1534/genetics.118.301468}.}
#'   \item{Gerard, D. and Ferrão, L. F. V. (2019). Priors for Genotyping Polyploids. \emph{Bioinformatics}. doi: \href{https://doi.org/10.1093/bioinformatics/btz852}{10.1093/bioinformatics/btz852}.}
#' }
#'
#' @export
#'
#' @author David Gerard
#'
#' @examples
#' ## F1 Population where parent 1 has 1 copy of the referenc allele
#' ## and parent 2 has 4 copies of the reference allele.
#' ploidy <- 6
#' rgeno(n = 10, ploidy = ploidy, model = "f1", p1geno = 1, p2geno = 4)
#'
#' ## A population in Hardy-Weinberge equilibrium with an
#' ## allele frequency of 0.75
#' rgeno(n = 10, ploidy = ploidy, model = "hw", allele_freq = 0.75)
#'
rgeno <- function(n,
                  ploidy,
                  model       = c("hw", "bb", "norm", "ash", "f1", "s1", "f1pp", "s1pp", "flex", "uniform"),
                  allele_freq = NULL,
                  od          = NULL,
                  p1geno      = NULL,
                  p2geno      = NULL,
                  mode        = NULL,
                  pivec       = NULL,
                  mu          = NULL,
                  sigma       = NULL,
                  p1_pair_weights = NULL,
                  p2_pair_weights = NULL) {
  ## Check input ----------------------------------------------------------
  model <- match.arg(model)
  assertthat::are_equal(length(ploidy), 1)
  assertthat::are_equal(length(n), 1)
  assertthat::assert_that(ploidy >= 0)
  if (model == "hw") {
    stopifnot(!is.null(allele_freq))
    stopifnot(is.null(p1geno))
    stopifnot(is.null(p2geno))
    stopifnot(is.null(mode))
    stopifnot(is.null(pivec))
    stopifnot(is.null(od))
    stopifnot(is.null(mu))
    stopifnot(is.null(sigma))
    stopifnot(is.null(p1_pair_weights))
    stopifnot(is.null(p2_pair_weights))
  } else if (model == "ash") {
    stopifnot(is.null(allele_freq))
    stopifnot(is.null(p1geno))
    stopifnot(is.null(p2geno))
    stopifnot(!is.null(mode))
    stopifnot(!is.null(pivec))
    stopifnot(is.null(od))
    stopifnot(is.null(mu))
    stopifnot(is.null(sigma))
    stopifnot(is.null(p1_pair_weights))
    stopifnot(is.null(p2_pair_weights))
  } else if (model == "f1") {
    stopifnot(is.null(allele_freq))
    stopifnot(!is.null(p1geno))
    stopifnot(!is.null(p2geno))
    stopifnot(is.null(mode))
    stopifnot(is.null(pivec))
    stopifnot(is.null(od))
    stopifnot(is.null(mu))
    stopifnot(is.null(sigma))
    stopifnot(is.null(p1_pair_weights))
    stopifnot(is.null(p2_pair_weights))
    if (ploidy %% 2 != 0) {
      stop("rgeno: ploidy must be even when model = 'f1'.")
    }
  } else if (model == "s1") {
    stopifnot(is.null(allele_freq))
    stopifnot(!is.null(p1geno))
    stopifnot(is.null(p2geno))
    stopifnot(is.null(mode))
    stopifnot(is.null(pivec))
    stopifnot(is.null(od))
    stopifnot(is.null(mu))
    stopifnot(is.null(sigma))
    stopifnot(is.null(p1_pair_weights))
    stopifnot(is.null(p2_pair_weights))
    if (ploidy %% 2 != 0) {
      stop("rgeno: ploidy must be even when model = 's1'.")
    }
  } else if (model == "f1pp") {
    stopifnot(is.null(allele_freq))
    stopifnot(!is.null(p1geno))
    stopifnot(!is.null(p2geno))
    stopifnot(is.null(mode))
    stopifnot(is.null(pivec))
    stopifnot(is.null(od))
    stopifnot(is.null(mu))
    stopifnot(is.null(sigma))
    stopifnot(!is.null(p1_pair_weights))
    stopifnot(!is.null(p2_pair_weights))
    if (ploidy %% 2 != 0) {
      stop("rgeno: ploidy must be even when model = 'f1'.")
    }
  } else if (model == "s1pp") {
    stopifnot(is.null(allele_freq))
    stopifnot(!is.null(p1geno))
    stopifnot(is.null(p2geno))
    stopifnot(is.null(mode))
    stopifnot(is.null(pivec))
    stopifnot(is.null(od))
    stopifnot(is.null(mu))
    stopifnot(is.null(sigma))
    stopifnot(!is.null(p1_pair_weights))
    stopifnot(is.null(p2_pair_weights))
    if (ploidy %% 2 != 0) {
      stop("rgeno: ploidy must be even when model = 's1'.")
    }
  } else if (model == "flex") {
    stopifnot(is.null(allele_freq))
    stopifnot(is.null(p1geno))
    stopifnot(is.null(p2geno))
    stopifnot(is.null(mode))
    stopifnot(!is.null(pivec))
    stopifnot(is.null(od))
    stopifnot(is.null(mu))
    stopifnot(is.null(sigma))
    stopifnot(is.null(p1_pair_weights))
    stopifnot(is.null(p2_pair_weights))
  } else if (model == "uniform") {
    stopifnot(is.null(allele_freq))
    stopifnot(is.null(p1geno))
    stopifnot(is.null(p2geno))
    stopifnot(is.null(mode))
    stopifnot(is.null(pivec))
    stopifnot(is.null(od))
    stopifnot(is.null(mu))
    stopifnot(is.null(sigma))
    stopifnot(is.null(p1_pair_weights))
    stopifnot(is.null(p2_pair_weights))
  } else if (model == "bb") {
    stopifnot(!is.null(allele_freq))
    stopifnot(is.null(p1geno))
    stopifnot(is.null(p2geno))
    stopifnot(is.null(mode))
    stopifnot(is.null(pivec))
    stopifnot(!is.null(od))
    stopifnot(is.null(mu))
    stopifnot(is.null(sigma))
    stopifnot(is.null(p1_pair_weights))
    stopifnot(is.null(p2_pair_weights))
  } else if (model == "norm") {
    stopifnot(is.null(allele_freq))
    stopifnot(is.null(p1geno))
    stopifnot(is.null(p2geno))
    stopifnot(is.null(mode))
    stopifnot(is.null(pivec))
    stopifnot(is.null(od))
    stopifnot(!is.null(mu))
    stopifnot(!is.null(sigma))
    stopifnot(is.null(p1_pair_weights))
    stopifnot(is.null(p2_pair_weights))
  } else{
    stop("rgeno: how did you get here?")
  }

  if (!is.null(allele_freq)) {
    stopifnot(length(allele_freq) == 1)
    stopifnot(allele_freq >= 0, allele_freq <= 1)
  }
  if (!is.null(od)) {
    stopifnot(length(od) == 1)
    stopifnot(od >= 0, od <= 1)
  }
  if (!is.null(p1geno)) {
    stopifnot(length(p1geno) == 1)
    stopifnot(p1geno >=0, p1geno <= ploidy)
    stopifnot(p1geno %% 1 == 0)
  }
  if (!is.null(p2geno)) {
    stopifnot(length(p2geno) == 1)
    stopifnot(p2geno >=0, p2geno <= ploidy)
    stopifnot(p2geno %% 1 == 0)
  }
  if (!is.null(mode)) {
    stopifnot(length(mode) == 1)
    if (mode %% 1 == 0) {
      stop("pgeno: mode cannot be an integer.")
    }
  }
  if (!is.null(pivec)) {
    stopifnot(length(pivec) == (ploidy + 1))
    stopifnot(abs(sum(pivec) - 1) < 10^-8)
    stopifnot(pivec >= 0, pivec <= 1)
  }
  if (!is.null(mu)) {
    stopifnot(length(mu) == 1)
  }
  if (!is.null(sigma)) {
    stopifnot(length(sigma) == 1)
    stopifnot(sigma > 0)
  }
  if (!is.null(p1_pair_weights)) {
    stopifnot(abs(sum(p1_pair_weights) - 1) < 10^-8)
    stopifnot(p1_pair_weights >= 0, p1_pair_weights <= 1)
  }
  if (!is.null(p2_pair_weights)) {
    stopifnot(abs(sum(p2_pair_weights) - 1) < 10^-8)
    stopifnot(p2_pair_weights >= 0, p2_pair_weights <= 1)
  }

  ## Generate probability vector ---------------------------------------

  if (model == "s1") {
    p2geno <- p1geno
  }
  if (model == "s1pp") {
    p2geno <- p1geno
    p2_pair_weights <- p1_pair_weights
  }

  if ((model == "s1") | (model == "f1")) {
    pivec <- get_q_array(ploidy = ploidy)[p1geno + 1, p2geno + 1, ]
  } else if (model == "hw") {
    pivec <- stats::dbinom(x = 0:ploidy, size = ploidy, prob = allele_freq)
  } else if (model == "bb") {
    pivec <- dbetabinom(x = 0:ploidy, size = ploidy, mu = allele_freq, rho = od, log = FALSE)
  } else if (model == "uniform") {
    pivec <- rep(x = 1 / (ploidy + 1), length = ploidy + 1)
  } else if (model == "flex") {
    pivec <- pivec
  } else if (model == "ash") {
    pivec <- get_probk_vec(pivec = pivec, model = "ash", mode = mode)
  } else if (model == "norm") {
    pivec <- stats::dnorm(x = 0:ploidy, mean = mu, sd = sigma, log = TRUE)
    pivec <- exp(pivec - log_sum_exp(pivec))
  } else if (model == "f1pp" | model == "s1pp") {
    blist <- get_bivalent_probs(ploidy = ploidy)
    which_p1 <- blist$lvec == p1geno
    which_p2 <- blist$lvec == p2geno
    if (sum(which_p1) != length(p1_pair_weights)) {
      stop("rgeno: p1_pair_weights should have length equal to the number\nof bivalent configurations in get_bivalent_probs() that corresponds\nto the dosage level of the parent (p1geno)\nas seen in lvec.")
    }
    if (sum(which_p2) != length(p2_pair_weights)) {
      stop("rgeno: p2_pair_weights should have length equal to the number\nof bivalent configurations in get_bivalent_probs() that corresponds\nto the dosage level of the parent (p2geno)\nas seen in lvec.")
    }
    p1segprob <- colSums(blist$probmat[which_p1, , drop = FALSE] * p1_pair_weights)
    p2segprob <- colSums(blist$probmat[which_p2, , drop = FALSE] * p2_pair_weights)
    pivec <- c(convolve_up(p1segprob, p2segprob))
  } else {
    stop("rgeno: how did you get here?")
  }

  ## Generate the genotypes -------------------------------------------
  genovec <- sample(x = 0:ploidy, size = n, replace = TRUE, prob = pivec)

  return(genovec)
}


#' Simulate GBS data from the \code{\link{flexdog}} likelihood.
#'
#' This will take a vector of genotypes and a vector of total read-counts, then generate a vector of reference
#' counts. To get the genotypes, you could use \code{\link{rgeno}}. The likelihood used to generate read-counts
#' is described in detail in Gerard et. al. (2018).
#'
#' @param sizevec A vector of total read-counts for the individuals.
#' @param geno A vector of genotypes for the individuals. I.e. the number of reference alleles
#'     each individual has.
#' @param ploidy The ploidy of the species.
#' @param seq The sequencing error rate.
#' @param bias The bias parameter. Pr(a read after selected) / Pr(A read after selected).
#' @param od The overdispersion parameter. See the
#'     Details of the \code{rho} variable in \code{\link{betabinom}}.
#'
#' @seealso \code{\link{rgeno}} for a way to generate genotypes of individuals. \code{\link{rbetabinom}}
#'     for how we generate the read-counts.
#'
#' @return A vector the same length as \code{sizevec}. The ith element
#'     is the number of reference counts for individual i.
#'
#' @references
#' \itemize{
#'   \item{Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018). Genotyping Polyploids from Messy Sequencing Data. \emph{Genetics}, 210(3), 789-807. doi: \href{https://doi.org/10.1534/genetics.118.301468}{10.1534/genetics.118.301468}.}
#'   \item{Gerard, D. and Ferrão, L. F. V. (2019). Priors for Genotyping Polyploids. \emph{Bioinformatics}. doi: \href{https://doi.org/10.1093/bioinformatics/btz852}{10.1093/bioinformatics/btz852}.}
#' }
#'
#' @export
#'
#' @author David Gerard
#'
#' @examples
#' set.seed(1)
#' n       <- 100
#' ploidy  <- 6
#'
#' ## Generate the genotypes of individuals from an F1 population,
#' ## where the first parent has 1 copy of the reference allele
#' ## and the second parent has two copies of the reference
#' ## allele.
#' genovec <- rgeno(n = n, ploidy = ploidy, model = "f1",
#'                  p1geno = 1, p2geno = 2)
#'
#' ## Get the total number of read-counts for each individual.
#' ## Ideally, you would take this from real data as the total
#' ## read-counts are definitely not Poisson.
#' sizevec <- stats::rpois(n = n, lambda = 200)
#'
#' ## Generate the counts of reads with the reference allele
#' ## when there is a strong bias for the reference allele
#' ## and there is no overdispersion.
#' refvec  <- rflexdog(sizevec = sizevec, geno = genovec,
#'                     ploidy = ploidy, seq = 0.001,
#'                     bias = 0.5, od = 0)
#'
#' ## Plot the simulated data using plot_geno.
#' plot_geno(refvec = refvec, sizevec = sizevec,
#'           ploidy = ploidy, seq = 0.001, bias = 0.5)
#'
rflexdog <- function(sizevec,
                     geno,
                     ploidy,
                     seq = 0.005,
                     bias = 1,
                     od = 0.001) {
  ## Check input -----------------------------------------------------------------------
  assertthat::are_equal(length(sizevec), length(geno))
  assertthat::are_equal(length(ploidy), length(seq), length(bias), length(od), 1)
  assertthat::assert_that(ploidy >= 0)
  assertthat::are_equal(ploidy %% 1, 0)
  stopifnot(geno <= ploidy)
  stopifnot(geno >= 0)
  stopifnot(sizevec >= 0)
  assertthat::assert_that(seq >= 0, seq <= 1)
  assertthat::assert_that(bias >= 0)
  assertthat::assert_that(od >= 0, od <= 1)

  ## Some boundaries for numerical stability
  bval <- 10 ^ -8
  if (od < bval) {
    od <- bval
  } else if (od > 1 - bval) {
    od <- 1 - bval
  }

  ## generate data
  xivec <- xi_fun(p = geno / ploidy, eps = seq, h = bias)
  refvec <- rbetabinom(n = length(sizevec), size = sizevec, mu = xivec, rho = rep(od, length = length(sizevec)))
  return(refvec)
}






