## Function to simulate from the flexdog model

#' Simulate individual genotypes from one of the supported \code{\link{flexdog}} models.
#'
#' This is simulate genotypes of a sample of individuals drawn from one of the populations supported by
#' the \code{\link{flexdog}}. Specificially, either generic populations, those with a unimodal genotype
#' distribution, those in Hardy-Weinberg equilibrium, those resulting from a biparental cross or a generation
#' of selfing, and populations with a uniform genotype distribution.
#'
#' The allowable variable values of \code{allele_freq}, \code{p1geno}, \code{p2geno}, \code{pivec},
#' and \code{mode} varies based on the value of \code{model}. If \code{model = "ash"}, then only
#' \code{mode} and \code{pivec} can be non-\code{NULL}. If \code{model = "flex"} then only
#' \code{pivec} can be non-\code{NULL}. If \code{model = "hw"}, then only \code{allele_freq} can
#' be non-\code{NULL}. If \code{model = "f1"} then only \code{p1geno} and \code{p2geno} can be non-\code{NULL}.
#' If \code{model = "s1"} then only \code{p1geno} can be non-\code{NULL}. If \code{model = "uniform"}, then
#' none of the above variables can be non-\code{NULL}.
#'
#' @param n The number of observations.
#' @param ploidy The ploidy of the species.
#' @param model What form should the prior take? Should the genotype
#'     distribution be unimodal (\code{"ash"}), generically
#'     any categorical distribution (\code{"flex"}), binomial as a
#'     result of assuming Hardy-Weinberg equilibrium (\code{"hw"}),
#'     a convolution of hypergeometics as a result that the population
#'     consists of either an F1 cross (\code{"f1"}) or an S1
#'     cross (\code{"s1"}), or fixed at a discrete uniform
#'     (\code{"uniform"})? See Details for more information.
#' @param allele_freq If \code{model = "hw"}, then this is the allele frequency of the population.
#'     For any other model, this should be \code{NULL}.
#' @param p1geno Either the first parent's genotype if \code{model = "f1"}, or the only parent's
#'     genotype if \code{model = "s1"}. For any other model, this should be \code{NULL}.
#' @param p2geno The second parent's genotype if \code{model = "f1"}.
#'     For any other model, this should be \code{NULL}.
#' @param pivec A vector of probabilities. If \code{model = "ash"}, then this represents
#'     the mixing proportions of the discrete uniforms. If
#'     \code{model = "flex"}, then element \code{i} is the probability of genotype \code{i - 1}.
#'     For any other model, this should be \code{NULL}.
#' @param mode If \code{model = "ash"}, this is the center of the distribution. This should be a non-integer value
#'     (so the mode is either the floor or the cieling of \code{mode}).
#'     For any other model, this should be \code{NULL}.
#'
#' @return A vector of length \code{n} with the genotypes of the sampled individuals.
#'
#' @export
#'
#' @author David Gerard
#'
rgeno <- function(n,
                  ploidy,
                  model       = c("hw", "ash", "f1", "s1", "flex", "uniform"),
                  allele_freq = NULL,
                  p1geno      = NULL,
                  p2geno      = NULL,
                  mode        = NULL,
                  pivec       = NULL) {
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
  } else if (model == "ash") {
    stopifnot(is.null(allele_freq))
    stopifnot(is.null(p1geno))
    stopifnot(is.null(p2geno))
    stopifnot(!is.null(mode))
    stopifnot(!is.null(pivec))
  } else if (model == "f1") {
    stopifnot(is.null(allele_freq))
    stopifnot(!is.null(p1geno))
    stopifnot(!is.null(p2geno))
    stopifnot(is.null(mode))
    stopifnot(is.null(pivec))
    if (ploidy %% 2 != 0) {
      stop("pgeno: ploidy must be even when model = 'f1'.")
    }
  } else if (model == "s1") {
    stopifnot(is.null(allele_freq))
    stopifnot(!is.null(p1geno))
    stopifnot(is.null(p2geno))
    stopifnot(is.null(mode))
    stopifnot(is.null(pivec))
    if (ploidy %% 2 != 0) {
      stop("pgeno: ploidy must be even when model = 's1'.")
    }
  } else if (model == "flex") {
    stopifnot(is.null(allele_freq))
    stopifnot(is.null(p1geno))
    stopifnot(is.null(p2geno))
    stopifnot(is.null(mode))
    stopifnot(!is.null(pivec))
  } else if (model == "uniform") {
    stopifnot(is.null(allele_freq))
    stopifnot(is.null(p1geno))
    stopifnot(is.null(p2geno))
    stopifnot(is.null(mode))
    stopifnot(is.null(pivec))
  } else{
    stop("rgeno: how did you get here?")
  }

  if (!is.null(allele_freq)) {
    stopifnot(length(allele_freq) == 1)
    stopifnot(allele_freq >= 0, allele_freq <= 1)
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
    stopifnot(sum(pivec) == 1)
    stopifnot(pivec >= 0, pivec <= 1)
  }

  ## Generate probability vector ---------------------------------------

  if (model == "s1") {
    p2geno <- p1geno
  }

  if ((model == "s1") | (model == "f1")) {
    pivec <- updog::get_q_array(ploidy = ploidy)[p1geno + 1, p2geno + 1, ]
  } else if (model == "hw") {
    pivec <- stats::dbinom(x = 0:ploidy, size = ploidy, prob = allele_freq)
  } else if (model == "uniform") {
    pivec <- rep(x = 1 / (ploidy + 1), length = ploidy + 1)
  } else if (model == "flex") {
    pivec <- pivec
  } else if (model == "ash") {
    pivec <- get_probk_vec(pivec = pivec, model = "ash", mode = mode)
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
#' counts. To get the genotypes, you could use \code{\link{rgeno}}.
#'
#' @param sizevec A vector of total read-counts for the individuals.
#' @param geno A vector of genotypes for the individuals. I.e. the number of reference alleles
#'     each individual has.
#' @param ploidy The ploidy of the species.
#' @param seq The sequencing error rate.
#' @param bias The bias parameter. Pr(a read after selected) / Pr(A read after selected).
#' @param od The overdispersion parameter. See Details of \code{rho} variable in \code{\link{betabinom}}.
#'
#' @seealso \code{\link{rgeno}} for a way to generate genotypes of individuals. \code{\link{rbetabinom}}
#'
#' @export
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
#' ## Plot the simulated data using plot_geno from the updog package.
#' updog::plot_geno(ocounts = refvec, osize = sizevec,
#'                  ploidy = ploidy, seq = 0.001, bias = 0.5)
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






