## Functions for ashdog and flexdog

#' Flexible genotyping for polyploids from next-generation sequencing data.
#'
#' Genotype polyploid individuals from next generation
#' sequencing (NGS) data while assuming the genotype distribution is one of
#' several forms. \code{flexdog} does this while accounting for allele bias,
#' overdispersion, sequencing error, and possibly outlying observations
#' (if \code{model = "f1"} or \code{model = "s1"}). The method is described
#' in detail in Gerard et. al. (2018) and Gerard and Ferrão (2019).
#'
#' @inherit flexdog_full
#'
#' @param bias_init A vector of initial values for the bias parameter
#'     over the multiple runs of \code{flexdog_full}.
#' @param ... Additional parameters to pass to \code{\link{flexdog_full}}.
#'
#' @seealso Run \code{browseVignettes(package = "updog")} in R for example usage.
#'     Other useful functions include:
#' \describe{
#'     \item{\code{\link{multidog}}}{For running \code{flexdog} on multiple
#'           SNP's.}
#'     \item{\code{\link{flexdog_full}}}{For additional parameter options
#'           when running \code{flexdog}.}
#'     \item{\code{\link{rgeno}}}{For simulating genotypes under the allowable
#'           prior models in \code{flexdog}.}
#'     \item{\code{\link{rflexdog}}}{For simulating read-counts under the
#'           assumed likelihood model in \code{flexdog}.}
#'     \item{\code{\link{plot.flexdog}}}{For plotting the output of
#'           \code{flexdog}.}
#'     \item{\code{\link{oracle_mis}}}{For calculating the oracle genotyping
#'           error rates. This is useful for read-depth calculations
#'           \emph{before} collecting data. After you have data, using
#'           the value of \code{prop_mis} is better.}
#'     \item{\code{\link{oracle_cor}}}{For calculating the correlation
#'           between the true genotypes and an oracle estimator
#'           (useful for read-depth calculations \emph{before}
#'           collecting data).}
#' }
#'
#' @author David Gerard
#'
#' @examples
#' \donttest{
#'
#' ## An S1 population where the first individual
#' ## is the parent. Fit assuming outliers.
#' data("snpdat")
#' ploidy  <- 6
#' refvec  <- snpdat$counts[snpdat$snp == "SNP3"]
#' sizevec <- snpdat$size[snpdat$snp == "SNP3"]
#' fout    <- flexdog(refvec   = refvec[-1],
#'                    sizevec  = sizevec[-1],
#'                    ploidy   = ploidy,
#'                    model    = "s1",
#'                    p1ref    = refvec[1],
#'                    p1size   = sizevec[1],
#'                    outliers = TRUE)
#' plot(fout)
#'
#' }
#'
#' ## A natural population. We will assume a
#' ## normal prior since there are so few
#' ## individuals.
#' data("uitdewilligen")
#' ploidy  <- 4
#' refvec  <- uitdewilligen$refmat[, 1]
#' sizevec <- uitdewilligen$sizemat[, 1]
#' fout    <- flexdog(refvec  = refvec,
#'                    sizevec = sizevec,
#'                    ploidy  = ploidy,
#'                    model   = "norm")
#' plot(fout)
#'
#'
#'
#' @export
flexdog <- function(refvec,
                    sizevec,
                    ploidy,
                    model       = c("norm", "hw", "bb", "ash", "s1", "s1pp", "f1", "f1pp", "flex", "uniform", "custom"),
                    p1ref       = NULL,
                    p1size      = NULL,
                    p2ref       = NULL,
                    p2size      = NULL,
                    snpname     = NULL,
                    bias_init   = exp(c(-1, -0.5, 0, 0.5, 1)),
                    verbose     = TRUE,
                    outliers    = FALSE,
                    prior_vec   = NULL,
                    ...) {
  assertthat::assert_that(all(bias_init > 0))
  model <- match.arg(model)

  if (verbose) {
    if ((length(refvec) < (10 * (ploidy + 1))) & (model == "flex")) {
      cat('Very few individuals for model = "flex"\nYou might want to try model = "norm" instead.\n\n')
    }
  }

  fout <- list()
  fout$llike <- -Inf
  for (bias_index in seq_along(bias_init)) {
    if (verbose) {
      cat("         Fit:", bias_index, "of", length(bias_init), "\n")
      cat("Initial Bias:", bias_init[bias_index], "\n")
    }

    fcurrent <- flexdog_full(refvec    = refvec,
                             sizevec   = sizevec,
                             ploidy    = ploidy,
                             model     = model,
                             p1ref     = p1ref,
                             p1size    = p1size,
                             p2ref     = p2ref,
                             p2size    = p2size,
                             snpname   = snpname,
                             bias      = bias_init[bias_index],
                             verbose   = FALSE,
                             outliers  = outliers,
                             prior_vec = prior_vec,
                             ...)

    if (verbose) {
      cat("Log-Likelihood:", fcurrent$llike, "\n")
    }

    if (fcurrent$llike > fout$llike) {
      fout <- fcurrent

      if (verbose) {
        cat("Keeping new fit.\n\n")
      }
    } else if (verbose) {
      cat("Keeping old fit.\n\n")
    }
  }

  if (verbose) {
    cat("Done!\n")
  }

  return(fout)
}

#' Flexible genotyping for polyploids from next-generation sequencing data.
#'
#' Genotype polyploid individuals from next generation
#' sequencing (NGS) data while assuming the genotype distribution is one of
#' several forms. \code{flexdog} does this while accounting for allele bias,
#' overdispersion, sequencing error, and possibly outlying observations
#' (if \code{model = "f1"} or \code{model = "s1"}). This function has more
#' options than \code{\link{flexdog}} and is only meant for expert users.
#' The method is described in detail in Gerard et. al. (2018) and
#' Gerard and Ferrão (2019).
#'
#' Possible values of the genotype distribution (values of \code{model}) are:
#' \describe{
#'   \item{\code{"norm"}}{A distribution whose genotype frequencies are proportional
#'       to the density value of a normal with some mean and some standard deviation.
#'       Unlike the \code{"bb"} and \code{"hw"} options, this will allow for
#'       distributions both more and less dispersed than a binomial.
#'       This seems to be the most robust to violations in modeling assumptions, and so is the
#'       default. This prior class was developed in Gerard and Ferrão (2019).}
#'   \item{\code{"hw"}}{A binomial distribution that results from assuming that
#'       the population is in Hardy-Weinberg equilibrium (HWE). This actually does
#'       pretty well even when there are minor to moderate deviations from HWE.
#'       Though it does not perform as well as the `"norm"` option when there
#'       are severe deviations from HWE.}
#'   \item{\code{"bb"}}{A beta-binomial distribution. This is an overdispersed
#'       version of \code{"hw"} and can be derived from a special case of the Balding-Nichols model.}
#'   \item{\code{"ash"}}{Any unimodal prior. This can sometimes be sensitive to violations
#'       in modeling assumptions, but tends to work better than the \code{"flex"} option.
#'       This prior class was developed in Gerard and Ferrão (2019).}
#'   \item{\code{"s1"}}{This prior assumes the individuals are
#'       all full-siblings resulting
#'       from one generation of selfing. I.e. there is only
#'       one parent. This model assumes
#'       a particular type of meiotic behavior: polysomic
#'       inheritance with
#'       bivalent, non-preferential pairing.
#'       Since this is a pretty strong and well-founded prior,
#'       we allow \code{outliers = TRUE} when \code{model = "s1"}.}
#'   \item{\code{"s1pp"}}{The same as \code{"s1"} but accounts for possible
#'       (and arbitrary levels of) preferential
#'       pairing during meiosis. The only supported values of \code{ploidy} right now for
#'       this option are \code{4} and \code{6}.}
#'   \item{\code{"f1"}}{This prior assumes the individuals are all
#'       full-siblings resulting
#'       from one generation of a bi-parental cross.
#'       This model assumes
#'       a particular type of meiotic behavior: polysomic
#'       inheritance with
#'       bivalent, non-preferential pairing. Since this is
#'       a pretty strong
#'       and well-founded prior, we allow \code{outliers = TRUE}
#'       when \code{model = "f1"}.}
#'   \item{\code{"f1pp"}}{The same as \code{"f1"} but accounts for possible
#'       (and arbitrary levels of) preferential
#'       pairing during meiosis. This option is mostly untested for values of \code{ploidy}
#'       greater than \code{6}.}
#'   \item{\code{"flex"}}{Generically any categorical distribution. Theoretically,
#'       this works well if you have a lot of individuals. In practice, it seems to
#'       be much less robust to violations in modeling assumptions.}
#'   \item{\code{"uniform"}}{A discrete uniform distribution. This should never
#'       be used in practice.}
#'   \item{\code{"custom"}}{A pre-specified prior distribution. You specify
#'       it using the \code{prior_vec} argument. You should almost never
#'       use this option in practice.}
#' }
#'
#' You might think a good default is \code{model = "uniform"} because it is
#' somehow an "uninformative prior." But it is very informative and tends to
#' work horribly in practice. The intuition is that it will estimate
#' the allele bias and sequencing error rates so that the estimated genotypes
#' are approximately uniform (since we are assuming that they are approximately
#' uniform). This will usually result in unintuitive genotyping since most
#' populations don't have a uniform genotype distribution.
#' I include it as an option only for completeness. Please don't use it.
#'
#' The value of \code{prop_mis} is a very intuitive measure for
#' the quality of the SNP. \code{prop_mis} is the posterior
#' proportion of individuals mis-genotyped. So if you want only SNPS
#' that accurately genotype, say, 95\% of the individuals, you could
#' discard all SNPs with a \code{prop_mis} over \code{0.05}.
#'
#' The value of \code{maxpostprob} is a very intuitive measure
#' for the quality of the genotype estimate of an individual.
#' This is the posterior probability of correctly genotyping
#' the individual when using \code{geno} (the posterior mode)
#' as the genotype estimate. So if you want to correctly genotype,
#' say, 95\% of individuals, you could discard all individuals
#' with a \code{maxpostprob} of under \code{0.95}. However, if you are
#' just going to impute missing genotypes later, you might consider
#' not discarding any individuals as \code{flexdog}'s genotype estimates will
#' probably be more accurate than other more naive approaches, such
#' as imputing using the grand mean.
#'
#' In most datasets I've examined, allelic bias is a major issue. However,
#' you may fit the model assuming no allelic bias by setting
#' \code{update_bias = FALSE} and \code{bias_init = 1}.
#'
#' Prior to using \code{flexdog}, during the read-mapping step,
#' you could try to get rid of allelic bias by
#' using WASP (\url{https://doi.org/10.1101/011221}). If you are successful
#' in removing the allelic bias (because its only source was the read-mapping
#' step), then setting \code{update_bias = FALSE} and \code{bias_init = 1}
#' would be reasonable. You can visually inspect SNPs for bias by
#' using \code{\link{plot_geno}}.
#'
#' \code{flexdog}, like most methods, is invariant to which allele you
#' label as the "reference" and which you label as the "alternative".
#' That is, if you set \code{refvec} with the number of alternative
#' read-counts, then the resulting genotype estimates
#' will be the estimated allele dosage of the alternative allele.
#'
#' @param refvec A vector of counts of reads of the reference allele.
#' @param sizevec A vector of total counts.
#' @param ploidy The ploidy of the species. Assumed to be the same for each
#'     individual.
#' @param model What form should the prior (genotype distribution) take?
#'     See Details for possible values.
#' @param verbose Should we output more (\code{TRUE}) or less
#'     (\code{FALSE})?
#' @param mean_bias The prior mean of the log-bias.
#' @param var_bias The prior variance of the log-bias.
#' @param mean_seq The prior mean of the logit of the sequencing
#'     error rate.
#' @param var_seq The prior variance of the logit of the sequencing
#'     error rate.
#' @param mean_od The prior mean of the logit of the overdispersion parameter.
#' @param var_od The prior variance of the logit of the overdispersion
#'     parameter.
#' @param seq The starting value of the sequencing error rate.
#' @param bias The starting value of the bias.
#' @param od The starting value of the overdispersion parameter.
#' @param mode The mode if \code{model = "ash"}. If not provided,
#'     \code{flexdog} will estimate the mode. This is the starting point of
#'     the allele frequency if \code{model = "hw"}. This should be
#'     \code{NULL} for all other options of \code{model}.
#' @param itermax The maximum number of EM iterations to run for each mode
#'     (if \code{model = "ash"}) or the total number of EM iterations to
#'     run (for any other value of \code{model}).
#' @param tol The tolerance stopping criterion. The EM algorithm will stop
#'     if the difference in the log-likelihoods between two consecutive
#'     iterations is less than \code{tol}.
#' @param update_bias A logical. Should we update \code{bias}
#'     (\code{TRUE}), or not (\code{FALSE})?
#' @param update_seq A logical. Should we update \code{seq}
#'     (\code{TRUE}), or not (\code{FALSE})?
#' @param update_od A logical. Should we update \code{od}
#'     (\code{TRUE}), or not (\code{FALSE})?
#' @param fs1_alpha The value at which to fix
#'     the mixing proportion for the uniform component
#'      when \code{model = "f1"}, \code{model = "f1pp"},
#'     \code{model = "s1"}, or \code{model = "s1pp"}.
#'     I would recommend some small
#'     value such as \code{10^-3}.
#' @param p1ref The reference counts for the first parent if
#'     \code{model = "f1"} (or \code{model = "f1pp"}), or for
#'     the only parent if \code{model = "s1"} (or \code{model = "s1pp"}).
#' @param p1size The total counts for the first parent if
#'     \code{model = "f1"} (or \code{model = "f1pp"}),
#'     or for the only parent if \code{model = "s1"} (or \code{model = "s1pp"}).
#' @param p2ref The reference counts for the second parent if
#'     \code{model = "f1"} (or \code{model = "f1pp"}).
#' @param p2size The total counts for the second parent if
#'     \code{model = "f1"} (or \code{model = "f1pp"}).
#' @param snpname A string. The name of the SNP under consideration.
#'     This is just returned in the \code{input} list for your reference.
#' @param ashpen The penalty to put on the unimodal prior.
#'     Larger values shrink the unimodal prior towards the
#'     discrete uniform distribution.
#' @param outliers A logical. Should we allow for the inclusion of outliers
#'     (\code{TRUE}) or not (\code{FALSE}). Only supported when
#'     \code{model = "f1"} or \code{model = "s1"}. I wouldn't
#'     recommend it for any other model anyway.
#' @param prior_vec The pre-specified genotype distribution. Only used if
#'     \code{model = "custom"} and must otherwise be \code{NULL}. If specified,
#'     then it should be a vector of length \code{ploidy + 1} with
#'     non-negative elements that sum to 1.
#'
#' @return An object of class \code{flexdog}, which consists
#'     of a list with some or all of the following elements:
#' \describe{
#'   \item{\code{bias}}{The estimated bias parameter.}
#'   \item{\code{seq}}{The estimated sequencing error rate.}
#'   \item{\code{od}}{The estimated overdispersion parameter.}
#'   \item{\code{num_iter}}{The number of EM iterations ran. You should
#'       be wary if this equals \code{itermax}.}
#'   \item{\code{llike}}{The maximum marginal log-likelihood.}
#'   \item{\code{postmat}}{A matrix of posterior probabilities of each
#'       genotype for each individual. The rows index the individuals
#'       and the columns index the allele dosage.}
#'   \item{\code{gene_dist}}{The estimated genotype distribution. The
#'       \code{i}th element is the proportion of individuals with
#'       genotype \code{i-1}. If \code{outliers = TRUE}, then this
#'       is conditional on the point not being an outlier.}
#'   \item{\code{par}}{A list of the final estimates of the parameters
#'       of the genotype distribution. The elements included in \code{par}
#'       depends on the value of \code{model}:
#'       \describe{
#'       \item{\code{model = "norm"}:}{\code{mu} is the normal mean and \code{sigma}
#'           is the normal standard deviation (not variance).}
#'       \item{\code{model = "hw"}:}{\code{alpha} is the major allele frequency.}
#'       \item{\code{model = "bb"}:}{\code{alpha} is the major allele frequency and
#'           \code{tau} is the overdispersion parameter (see the description of
#'           \code{rho} in the Details of \code{\link{betabinom}}).}
#'       \item{\code{model = "ash"}:}{\code{par} is an empty list.}
#'       \item{\code{model = "s1"}:}{\code{pgeno} is the allele dosage of the parent and
#'           \code{alpha} is the mixture proportion of the discrete uniform (included and
#'           fixed at a small value mostly for numerical stability reasons). See the description
#'           of \code{fs1_alpha} in \code{\link{flexdog_full}}.}
#'       \item{\code{model = "s1pp"}:}{\code{pgeno} is the allele dosage of the parent and
#'           \code{p1_pair_weights} contains a vector of mixing weights where element \code{i} is the
#'           mixing proportion for the segregation distribution in row \code{i} of
#'           \code{get_bivalent_probs(ploidy)$probmat[get_bivalent_probs(ploidy)$lvec == pgeno, , drop = FALSE]}.}
#'       \item{\code{model = "f1"}:}{\code{p1geno} is the allele dosage of the first parent,
#'           \code{p2geno} is the allele dosage of the second parent, and \code{alpha}
#'           is the mixture proportion of the discrete uniform (included and
#'           fixed at a small value mostly for numerical stability reasons). See the description
#'           of \code{fs1_alpha} in \code{\link{flexdog_full}}.}
#'       \item{\code{model = "f1pp"}:}{\code{p1geno} is the allele dosage of the first parent,
#'           \code{p2geno} is the allele dosage of the second parent,
#'           \code{p1_pair_weights} contains a vector of mixing weights where element \code{i} is the
#'           mixing proportion for the segregation distribution for parent 1 in row \code{i} of
#'           \code{get_bivalent_probs(ploidy)$probmat[get_bivalent_probs(ploidy)$lvec == p1geno, , drop = FALSE]},
#'           and \code{p2_pair_weights} contains a vector of mixing weights where element \code{i} is the
#'           mixing proportion for the segregation distribution for parent 2 in row \code{i} of
#'           \code{get_bivalent_probs(ploidy)$probmat[get_bivalent_probs(ploidy)$lvec == p2geno, , drop = FALSE]}.}
#'       \item{\code{model = "flex"}:}{\code{par} is an empty list.}
#'       \item{\code{model = "uniform"}:}{\code{par} is an empty list.}
#'       \item{\code{model = "custom"}:}{\code{par} is an empty list.}
#'       }}
#'   \item{\code{geno}}{The posterior mode genotype. These are your
#'       genotype estimates.}
#'   \item{\code{maxpostprob}}{The maximum posterior probability. This
#'       is equivalent to the posterior probability of correctly
#'       genotyping each individual.}
#'   \item{\code{postmean}}{The posterior mean genotype. In downstream
#'       association studies, you might want to consider using these
#'       estimates.}
#'   \item{\code{input$refvec}}{The value of \code{refvec} provided by
#'       the user.}
#'   \item{\code{input$sizevec}}{The value of \code{sizevec} provided
#'       by the user.}
#'   \item{\code{input$ploidy}}{The value of \code{ploidy} provided
#'       by the user.}
#'   \item{\code{input$model}}{The value of \code{model} provided by
#'       the user.}
#'   \item{\code{input$p1ref}}{The value of \code{p1ref} provided by the user.}
#'   \item{\code{input$p1size}}{The value of \code{p1size} provided by the user.}
#'   \item{\code{input$p2ref}}{The value of \code{p2ref} provided by the user.}
#'   \item{\code{input$p2size}}{The value of \code{p2size} provided by the user.}
#'   \item{\code{input$snpname}}{The value of \code{snpname} provided by the user.}
#'   \item{\code{prop_mis}}{The posterior proportion of individuals
#'       genotyped incorrectly.}
#'   \item{\code{out_prop}}{The estimated proportion of points that
#'       are outliers. Only available if \code{outliers = TRUE}.}
#'   \item{\code{prob_out}}{The ith element is the posterior probability
#'       that individual i is an outlier. Only available if
#'       \code{outliers = TRUE}.}
#' }
#'
#' @author David Gerard
#'
#' @references
#' \itemize{
#'   \item{Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018). Genotyping Polyploids from Messy Sequencing Data. \emph{Genetics}, 210(3), 789-807. doi: \href{https://doi.org/10.1534/genetics.118.301468}{10.1534/genetics.118.301468}.}
#'   \item{Gerard, D. and Ferrão, L. F. V. (2019). Priors for Genotyping Polyploids. \emph{Bioinformatics}. doi: \href{https://doi.org/10.1093/bioinformatics/btz852}{10.1093/bioinformatics/btz852}.}
#' }
#'
#' @seealso
#' Run \code{browseVignettes(package = "updog")} in R for example usage.
#'     Other useful functions include:
#' \describe{
#'     \item{\code{\link{multidog}}}{For running \code{flexdog} on multiple
#'           SNP's.}
#'     \item{\code{\link{flexdog}}}{For a more user-friendly version of
#'           \code{flexdog_full}.}
#'     \item{\code{\link{rgeno}}}{For simulating genotypes under the allowable
#'           prior models in \code{flexdog}.}
#'     \item{\code{\link{rflexdog}}}{For simulating read-counts under the
#'           assumed likelihood model in \code{flexdog}.}
#'     \item{\code{\link{plot.flexdog}}}{For plotting the output of
#'           \code{flexdog}.}
#'     \item{\code{\link{oracle_mis}}}{For calculating the oracle genotyping
#'           error rates. This is useful for read-depth calculations
#'           \emph{before} collecting data. After you have data, using
#'           the value of \code{prop_mis} is better.}
#'     \item{\code{\link{oracle_cor}}}{For calculating the correlation
#'           between the true genotypes and an oracle estimator
#'           (useful for read-depth calculations \emph{before}
#'           collecting data).}
#' }
#'
#' @examples
#' ## A natural population. We will assume a
#' ## normal prior since there are so few
#' ## individuals.
#' data("uitdewilligen")
#' ploidy  <- 4
#' refvec  <- uitdewilligen$refmat[, 1]
#' sizevec <- uitdewilligen$sizemat[, 1]
#' fout    <- flexdog_full(refvec  = refvec,
#'                         sizevec = sizevec,
#'                         ploidy  = ploidy,
#'                         model   = "norm")
#' plot(fout)
#'
#' @export
flexdog_full <- function(refvec,
                         sizevec,
                         ploidy,
                         model       = c("norm", "hw", "bb", "ash", "s1", "s1pp", "f1", "f1pp", "flex", "uniform", "custom"),
                         verbose     = TRUE,
                         mean_bias   = 0,
                         var_bias    = 0.7 ^ 2,
                         mean_seq    = -4.7,
                         var_seq     = 1,
                         mean_od     = -5.5,
                         var_od      = 0.5 ^ 2,
                         seq         = 0.005,
                         bias        = 1,
                         od          = 0.001,
                         update_bias = TRUE,
                         update_seq  = TRUE,
                         update_od   = TRUE,
                         mode        = NULL,
                         itermax     = 200,
                         tol         = 10 ^ -4,
                         fs1_alpha   = 10 ^ -3,
                         ashpen      = 10 ^ -6,
                         p1ref       = NULL,
                         p1size      = NULL,
                         p2ref       = NULL,
                         p2size      = NULL,
                         snpname     = NULL,
                         outliers    = FALSE,
                         prior_vec   = NULL) {

  ## Check input -----------------------------------------------------
  model <- match.arg(model)
  if (model == "uniform") {
    warning("flexdog: Using model = 'uniform' is almost always a bad idea.\nTry model = 'hw' or model = 'norm' if you have data from a population study.")
  }
  if (outliers) {
    if ((model != "s1") & (model != "f1")) {
      stop('flexdog: outliers = TRUE only supported when model = "f1" or model = "s1".')
    }
  }
  if ((model == "f1pp" | model == "s1pp" | model == "f1ppdr" | model == "s1ppdr") & ploidy == 2) {
    stop("flexdog: preferential pairing and double reduction\ncannot occur in diploids (i.e. when ploidy = 2).\nTry using f1 or s1 instead.")
  }
  if ((model == "s1pp") & (ploidy != 6) & (ploidy != 4)) {
    stop("flexdog: model = 's1pp' is only supported for ploides 4 and 6.")
  }
  if (model == "s1ppdr") {
    stop("flexdog: model = 's1ppdr' is not supported right now.")
  }
  if (model == "f1ppdr") {
    warning("flexdog: model = 'f1ppdr' is a mostly untested option.\nUse at your own risk.")
  }
  if (model == "f1pp" & ploidy > 6) {
    warning("flexdog: flexdog has not been extensively tested\nfor model = 'f1pp' when ploidy > 6.\nUse at your own risk.")
  }

  assertthat::are_equal(length(refvec), length(sizevec))
  assertthat::assert_that(all(sizevec >= refvec, na.rm = TRUE))
  assertthat::assert_that(all(refvec >= 0, na.rm = TRUE))
  assertthat::are_equal(1, length(verbose), length(mean_bias),
                        length(var_bias), length(mean_seq),
                        length(var_seq), length(seq),
                        length(bias), length(od),
                        length(mean_od), length(var_od))
  assertthat::assert_that(var_bias > 0)
  assertthat::assert_that(var_seq > 0)
  assertthat::assert_that(seq >= 0, seq <= 1)
  assertthat::assert_that(bias > 0)
  assertthat::assert_that(od >= 0, od <= 1)
  assertthat::assert_that(is.logical(verbose))
  assertthat::are_equal(ploidy %% 1, 0)
  assertthat::assert_that(ploidy > 0)
  assertthat::assert_that(tol > 0)
  assertthat::are_equal(itermax %% 1, 0)
  assertthat::assert_that(itermax > 0)
  assertthat::assert_that(is.logical(update_bias))
  assertthat::assert_that(is.logical(update_seq))
  assertthat::assert_that(is.logical(update_od))
  assertthat::assert_that(ashpen >= 0)

  ## check fs1_alpha -----------------------------------------------
  assertthat::are_equal(length(fs1_alpha), 1)
  assertthat::assert_that(fs1_alpha <= 1, fs1_alpha >= 0)

  ## Check custom --------------------------------------------------
  if (!is.null(prior_vec) & model != "custom") {
    stop('Cannot specify `prior_vec` when model != "custom"')
  } else if (!is.null(prior_vec)) {
    stopifnot(prior_vec >= 0)
    stopifnot(prior_vec <= 1)
    stopifnot(abs(sum(prior_vec) - 1) < 10^-6)
    stopifnot(length(prior_vec) == ploidy + 1)
  } else if (model == "custom" & is.null(prior_vec)) {
    stop('must specify `prior_vec` when model = "custom"')
  }

  ## Check p1ref, p2ref, p1size, p2size ----------------------------
  if ((!is.null(p1ref) | !is.null(p1size)) & (model != "f1" & model != "s1" & model != "f1pp" & model != "s1pp" & model != "f1ppdr" & model != "s1ppdr")) {
    stop("flexdog: if model is not 'f1', 's1', 'f1pp', 's1pp', 'f1ppdr', or 's1ppdr' then p1ref and p1size both need to be NULL.")
  }
  if ((!is.null(p2ref) | !is.null(p2size)) & (model != "f1" & model != "f1pp" & model != "f1ppdr")) {
    stop("flexdog: if model is not 'f1', 'f1pp', or 'f1ppdr' then p2ref and p2size both need to be NULL.")
  }
  if ((is.null(p1ref) & !is.null(p1size)) | (!is.null(p1ref) & is.null(p1size))) {
    stop("flexdog: p1ref and p1size either need to be both NULL or both non-NULL.")
  }
  if ((is.null(p2ref) & !is.null(p2size)) | (!is.null(p2ref) & is.null(p2size))) {
    stop("flexdog: p1ref and p1size either need to be both NULL or both non-NULL.")
  }
  if (!is.null(p1ref)) {
    stopifnot(p1ref >= 0, p1size >= p1ref)
  }
  if (!is.null(p2ref)) {
    stopifnot(p2ref >= 0, p2size >= p2ref)
  }
  if (!is.null(snpname)) {
    stopifnot(is.character(snpname))
    stopifnot(length(snpname) == 1)
  }

  ## check and set mode under various models -----------------------
  if (!is.null(mode) & model == "flex") {
    stop('flexdog: `model` cannot equal `"flex"` when `mode` is specified.')
  } else if (is.null(mode) & model == "flex") {
    mode_vec <- 0
  } else if (!is.null(mode) & (model == "f1" | model == "s1" | model == "uniform" | model == "f1pp" | model == "s1pp" | model == "f1ppdr" | model == "s1ppdr" | model == "custom")) {
    stop('flexdog: `model` cannot equal "f1", "s1", "f1pp", "s1pp", "f1ppdr", "s1ppdr", "uniform", or "custom" when `mode` is specified.')
  } else if (is.null(mode) & (model == "f1" | model == "s1" | model == "uniform" | model == "f1pp" | model == "s1pp" | model == "f1ppdr" | model == "s1ppdr")) {
    mode_vec <- mean(refvec / sizevec, na.rm = TRUE) ## just to initialize pivec
  } else if (!is.null(mode) & model == "ash") {
    stopifnot(length(mode) == 1)
    mode_vec <- mode
  } else if (is.null(mode) & model == "ash") {
    mode_vec <- (0:(ploidy - 1)) + 0.5
  } else if (is.null(mode) & (model == "hw" | model == "bb" | model == "norm")) {
    mode_vec <- mean(refvec / sizevec, na.rm = TRUE)
  } else if (!is.null(mode) & (model == "hw" | model == "bb" | model == "norm")) {
    if (any((mode < 0) | (mode > 1))) {
      stop('If model = "hw" or model = "bb" or model = "norm" then `mode` should be between 0 and 1.\nIt is the initialization of the allele frequency.')
    }
  } else if (model == "custom") {
    mode_vec <- sum((0:ploidy) * prior_vec) / ploidy
  } else {
    stop("flexdog: Checking mode. How did you get here?")
  }

  ## Deal with missingness in sizevec and refvec -----------------------
  not_na_vec  <- !(is.na(refvec) | is.na(sizevec))
  refvec      <- refvec[not_na_vec]
  sizevec     <- sizevec[not_na_vec]

  ## Some variables needed to run EM ---------------------------
  control <- list() ## will contain parameters used to update pivec
  boundary_tol <- 10 ^ -6 ## smallest values of seq, bias, od,
                          ## how close can od and seq get to 1.
  if (od < boundary_tol) {
    od <- boundary_tol
  }
  if (seq < boundary_tol) {
    seq <- boundary_tol
  }
  if (bias < boundary_tol) {
    bias <- boundary_tol
  }
  if (od > 1 - boundary_tol) {
    od <- 1 - boundary_tol
  }
  if (seq > 1 - boundary_tol) {
    seq <- 1 - boundary_tol
  }
  if (ashpen < boundary_tol) {
    ashpen <- boundary_tol
  }

  ## since re-write these
  bias_init <- bias
  seq_init  <- seq
  od_init   <- od

  ## Run EM for each mode in `mode_vec` -----------------------
  return_list <- list(llike = -Inf)
  for (em_index in seq_along(mode_vec)) {
    mode <- mode_vec[em_index]
    bias <- bias_init
    seq  <- seq_init
    od   <- od_init

    if (verbose) {
      cat("Mode:", mode, "\n")
    }

    ## Get inner weight vec only once
    ## Used in convex optimization program
    if (model == "ash") {
      control$inner_weights <- get_inner_weights(ploidy = ploidy, mode = mode)
      control$lambda        <- ashpen
    } else if (model == "f1" | model == "s1") {
      control$qarray    <- get_q_array(ploidy = ploidy)
      control$fs1_alpha <- fs1_alpha
      control$outliers  <- outliers
      if (outliers) {
        control$out_prop  <- 0.01
      }
    } else if (model == "bb") {
      control$alpha <- mode ## initialize allele frequency for bb
      control$tau   <- boundary_tol ## initialize od for bb
    } else if (model == "norm") {
      control$mu    <- mode * ploidy ## initialize mean of normal
      control$sigma <- sqrt(ploidy * mode * (1 - mode)) ## initialize sd of normal
    } else if (model == "f1pp" | model == "s1pp") {
      control$blist <- get_bivalent_probs(ploidy = ploidy)
      control$fs1_alpha <- fs1_alpha
      control$p1_pair_weights <- list()
      for (ell in 0:ploidy) { ## initialize bivalent pairing weights
        control$p1_pair_weights[[ell + 1]] <- get_hyper_weights(ploidy = ploidy, ell = ell)$weightvec
      }
      if (model == "f1pp") {
        control$p2_pair_weights <- list()
        for (ell in 0:ploidy) {
          control$p2_pair_weights[[ell + 1]] <- get_hyper_weights(ploidy = ploidy, ell = ell)$weightvec
        }
      }
    } else if (model == "f1ppdr" | model == "s1ppdr") {
      control$blist <- get_bivalent_probs_dr(ploidy = ploidy)
      control$fs1_alpha <- fs1_alpha
      control$p1_pair_weights <- list()
      for (ell in 0:ploidy) { ## initialize bivalent pairing weights
        control$p1_pair_weights[[ell + 1]] <- rep(0, length = sum(control$blist$lvec == ell))
        control$p1_pair_weights[[ell + 1]][control$blist$penvec[control$blist$lvec == ell]] <- get_hyper_weights(ploidy = ploidy, ell = ell)$weightvec
        control$p1_pair_weights[[ell + 1]] <- control$p1_pair_weights[[ell + 1]] * 0.99 + 0.01 / length(control$p1_pair_weights[[ell + 1]])
      }
      if (model == "f1ppdr") {
        control$p2_pair_weights <- list()
        for (ell in 0:ploidy) {
          control$p2_pair_weights[[ell + 1]] <- rep(0, length = sum(control$blist$lvec == ell))
          control$p2_pair_weights[[ell + 1]][control$blist$penvec[control$blist$lvec == ell]] <- get_hyper_weights(ploidy = ploidy, ell = ell)$weightvec
          control$p2_pair_weights[[ell + 1]] <- control$p2_pair_weights[[ell + 1]] * 0.99 + 0.01 / length(control$p2_pair_weights[[ell + 1]])
        }
      }
      control$mixing_pen <- control$blist$penvec * 1 ## The penalties to use on the mixing weights of the ppdr model.
    } else if (model == "custom") {
      control$prior_vec <- prior_vec
    }

    ## Initialize pivec so that two modes have equal prob if model = "ash".
    ##     Uniform if model = "flex".
    pivec <- initialize_pivec(ploidy = ploidy, mode = mode, model = model)
    if (model == "custom") { ## hack to get around passing more arguments to initialize_pivec
      pivec <- prior_vec
    }
    assertthat::are_equal(sum(pivec), 1)
    control$pivec <- pivec ## initialization for unimodal optimization

    probk_vec <- get_probk_vec(pivec = pivec, model = model, mode = mode)
    assertthat::are_equal(sum(probk_vec), 1)

    ## Run EM ----------------------------------------
    iter_index  <- 1
    err         <- tol + 1
    llike       <- -Inf
    while (err > tol & iter_index <= itermax) {
      llike_old <- llike

      ## E-step ----------------------
      if (!outliers) {
        wik_mat <- get_wik_mat(probk_vec = probk_vec,
                               refvec    = refvec,
                               sizevec   = sizevec,
                               ploidy    = ploidy,
                               seq       = seq,
                               bias      = bias,
                               od = od)
      } else {
        wik_temp <- get_wik_mat_out(probk_vec = probk_vec,
                                    out_prop  = control$out_prop,
                                    refvec    = refvec,
                                    sizevec   = sizevec,
                                    ploidy    = ploidy,
                                    seq       = seq,
                                    bias      = bias,
                                    od        = od)
        wik_mat <- wik_temp[, seq_len(ploidy + 1), drop = FALSE]
        prob_outlier <- wik_temp[, ploidy + 2]
      }

      ## Update seq, bias, and od ----
      oout <- stats::optim(par         = c(seq, bias, od),
                           fn          = obj_for_eps,
                           gr          = grad_for_eps,
                           method      = "L-BFGS-B",
                           lower       = rep(boundary_tol, 3),
                           upper       = c(1 - boundary_tol, Inf,
                                           1 - boundary_tol),
                           control     = list(fnscale = -1, maxit = 20),
                           refvec      = refvec,
                           sizevec     = sizevec,
                           ploidy      = ploidy,
                           mean_bias   = mean_bias,
                           var_bias    = var_bias,
                           mean_seq    = mean_seq,
                           var_seq     = var_seq,
                           mean_od     = mean_od,
                           var_od      = var_od,
                           wmat        = wik_mat,
                           update_seq  = update_seq,
                           update_bias = update_bias,
                           update_od   = update_od)
      seq  <- oout$par[1]
      bias <- oout$par[2]
      od   <- oout$par[3]

      ## if F1 or S1, update betabinomial log-likelihood of parent counts
      if (model == "f1" | model == "f1pp" | model == "f1ppdr") {
        if (!is.null(p1ref)) {
          xi_vec <- xi_fun(p = 0:ploidy / ploidy, eps = seq, h = bias)
          control$p1_lbb <- dbetabinom(x    = rep(p1ref, ploidy + 1),
                                       size = rep(p1size, ploidy + 1),
                                       mu   = xi_vec,
                                       rho  = od,
                                       log  = TRUE)
        }
        if (!is.null(p2ref)) {
          xi_vec <- xi_fun(p = 0:ploidy / ploidy, eps = seq, h = bias)
          control$p2_lbb <- dbetabinom(x    = rep(p2ref, ploidy + 1),
                                       size = rep(p2size, ploidy + 1),
                                       mu   = xi_vec,
                                       rho  = od,
                                       log  = TRUE)
        }
      } else if (model == "s1" | model == "s1pp" | model == "s1ppdr") {
        if (!is.null(p1ref)) {
          xi_vec <- xi_fun(p = 0:ploidy / ploidy, eps = seq, h = bias)
          control$p1_lbb <- dbetabinom(x    = rep(p1ref, ploidy + 1),
                                       size = rep(p1size, ploidy + 1),
                                       mu   = xi_vec,
                                       rho  = od,
                                       log  = TRUE)
        }
      } else {
        ## do nothing
      }

      ## Update pivec ----------------
      weight_vec <- colSums(wik_mat)
      if (outliers) {
        control$weight_out <- sum(prob_outlier)
      }
      fupdate_out <- flex_update_pivec(weight_vec = weight_vec,
                                       model      = model,
                                       control    = control)
      pivec <- fupdate_out$pivec
      control$pivec <- pivec ## initial condition for unimodal optimization

      ## New initialization parameters for priors that use gradient ascent.
      if (model == "bb") {
        control$alpha <- fupdate_out$par$alpha
        control$tau   <- fupdate_out$par$tau
      } else if (model == "norm") {
        control$mu    <- fupdate_out$par$mu
        control$sigma <- fupdate_out$par$sigma
      } else if (((model == "f1") | (model == "s1")) & outliers) {
        out_prop         <- fupdate_out$par$out_prop
        control$out_prop <- fupdate_out$par$out_prop
      } else if (model == "f1pp" | model == "f1ppdr") {
        control$p1_pair_weights[[fupdate_out$par$p1geno + 1]] <-
          fupdate_out$par$p1_pair_weights
        control$p2_pair_weights[[fupdate_out$par$p2geno + 1]] <-
          fupdate_out$par$p2_pair_weights
      } else if (model == "s1pp" | model == "s1ppdr") {
        control$p1_pair_weights[[fupdate_out$par$p1geno + 1]] <-
          fupdate_out$par$p1_pair_weights
      }

      ## Update probk_vec -----------------------------------------------
      pivec[pivec < 0] <- 0
      pivec[pivec > 1] <- 1
      probk_vec <- get_probk_vec(pivec = pivec, model = model, mode = mode)

      ## Calculate likelihood and update stopping criteria --------------
      if (!outliers) {
        llike <- flexdog_obj(probk_vec = probk_vec,
                             refvec    = refvec,
                             sizevec   = sizevec,
                             ploidy    = ploidy,
                             seq       = seq,
                             bias      = bias,
                             od        = od,
                             mean_bias = mean_bias,
                             var_bias  = var_bias,
                             mean_seq  = mean_seq,
                             var_seq   = var_seq,
                             mean_od   = mean_od,
                             var_od    = var_od)
      } else {
        llike <- flexdog_obj_out(probk_vec = probk_vec,
                                 out_prop  = out_prop,
                                 refvec    = refvec,
                                 sizevec   = sizevec,
                                 ploidy    = ploidy,
                                 seq       = seq,
                                 bias      = bias,
                                 od        = od,
                                 mean_bias = mean_bias,
                                 var_bias  = var_bias,
                                 mean_seq  = mean_seq,
                                 var_seq   = var_seq,
                                 mean_od   = mean_od,
                                 var_od    = var_od)
      }


      if (model == "ash") { ## add small penalty if "ash"
        llike <- llike + ashpen_fun(lambda = ashpen, pivec = pivec)
      } else if (model == "f1ppdr") { ## add more penalties if f1ppdr
        llike <- llike +
          dr_pen(pairweights = fupdate_out$par$p1_pair_weights,
                 mixing_pen = control$mixing_pen[control$blist$lvec == fupdate_out$par$p1geno]) +
          dr_pen(pairweights = fupdate_out$par$p2_pair_weights,
                 mixing_pen = control$mixing_pen[control$blist$lvec == fupdate_out$par$p2geno])
      } else if (model == "s1ppdr") { ## NEED TO FIX if I ever finally support s1ppdr
        llike <- llike +
          dr_pen(pairweights = fupdate_out$par$p1_pair_weights,
                 mixing_pen = control$mixing_pen[control$blist$lvec == fupdate_out$par$p1geno])
      }

      err        <- abs(llike - llike_old)
      iter_index <- iter_index + 1

      if (llike < llike_old - 10 ^ -5) {
        warning(paste0("flexdog: likelihood not increasing.\nDifference is: ",
                       llike - llike_old))
        cat("Warn:\n")
        cat("Diff:", llike - llike_old, "\n")
        cat("p1geno:", fupdate_out$par$p1geno, "\n")
        cat("p2geno:", fupdate_out$par$p2geno, "\n")
        cat("llike:", format(llike, digits = 10), "\n")
        cat("llike_old:", format(llike_old, digits = 10), "\n\n")
        cat("pivec:", pivec, "\n")
        if (verbose) {
          cat("\nindex: ", iter_index, "\n")
          cat("llike: ", llike, "\n")
          cat("pivec: ", paste0("c(", paste0(pivec, collapse = ", "), ")"), "\n")
          cat("weight_vec: ",  paste0("c(", paste0(weight_vec, collapse = ", "), ")"), "\n\n")
        }
      }
    }

    if (verbose) {
      cat("llike:", llike, "\n\n")
    }

    ## Check which mode has the highest likelihood --------------------
    temp_list <- list(bias      = bias,
                      seq       = seq,
                      od        = od,
                      num_iter  = iter_index,
                      llike     = llike,
                      postmat   = wik_mat,
                      gene_dist = probk_vec,
                      par       = fupdate_out$par)
    if (outliers) {
      temp_list$prob_outlier    <- prob_outlier
      temp_list$out_prop        <- out_prop
    }
    if (temp_list$llike > return_list$llike) {
      return_list <- temp_list
    }
  }

  ## Adjust postmat for outliers --------------------------------------
  if (outliers) {
    return_list$postmat <- return_list$postmat +
      outer(return_list$prob_outlier, return_list$gene_dist, "*")

    temp                     <- rep(NA, length = length(not_na_vec))
    temp[not_na_vec]         <- return_list$prob_outlier
    return_list$prob_outlier <- temp
  }


  ## Summaries --------------------------------------------------------
  return_list$geno          <- apply(return_list$postmat, 1, which.max) - 1
  return_list$maxpostprob   <- return_list$postmat[cbind(seq_len(nrow(return_list$postmat)), return_list$geno + 1)]
  return_list$postmean      <- c(return_list$postmat %*% 0:ploidy)
  return_list$input$refvec  <- refvec
  return_list$input$sizevec <- sizevec
  return_list$input$ploidy  <- ploidy
  return_list$input$model   <- model
  return_list$input$p1ref   <- p1ref
  return_list$input$p1size  <- p1size
  return_list$input$p2ref   <- p2ref
  return_list$input$p2size  <- p2size
  return_list$input$snpname <- snpname
  return_list$prop_mis      <- 1 - mean(return_list$maxpostprob)

  ## Add back missingness ---------------------------------------------
  temp                     <- rep(NA, length = length(not_na_vec))
  temp[not_na_vec]         <- refvec
  return_list$input$refvec <- temp

  temp                      <- rep(NA, length = length(not_na_vec))
  temp[not_na_vec]          <- sizevec
  return_list$input$sizevec <- temp

  temp             <- rep(NA, nrow = length(not_na_vec))
  temp[not_na_vec] <- return_list$geno
  return_list$geno <- temp

  temp                    <- rep(NA, nrow = length(not_na_vec))
  temp[not_na_vec]        <- return_list$maxpostprob
  return_list$maxpostprob <- temp

  temp                 <- rep(NA, nrow = length(not_na_vec))
  temp[not_na_vec]     <- return_list$postmean
  return_list$postmean <- temp

  temp                <- matrix(NA, nrow = length(not_na_vec),
                                ncol = ncol(return_list$postmat))
  temp[not_na_vec, ]  <- return_list$postmat
  return_list$postmat <- temp

  ## Set class to flexdog ---------------------------------------------
  class(return_list) <- "flexdog"

  return(return_list)
}

#' Draw a genotype plot from the output of \code{\link{flexdog}}.
#'
#' @inherit plot.mupdog description details
#'
#' @param x A \code{flexdog} object.
#' @param use_colorblind Should we use a colorblind-safe palette
#'     (\code{TRUE}) or not (\code{FALSE})? \code{TRUE} is only allowed
#'     if the ploidy is less than or equal to 6.
#' @param ... Not used.
#'
#' @author David Gerard
#'
#' @seealso
#' \describe{
#'   \item{\code{\link{plot_geno}}}{The underlying plotting function.}
#'   \item{\code{\link{flexdog}}}{Creates a \code{flexdog} object.}
#' }
#'
#' @references
#' \itemize{
#'   \item{Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018). Genotyping Polyploids from Messy Sequencing Data. \emph{Genetics}, 210(3), 789-807. doi: \href{https://doi.org/10.1534/genetics.118.301468}{10.1534/genetics.118.301468}.}
#'   \item{Gerard, D. and Ferrão, L. F. V. (2019). Priors for Genotyping Polyploids. \emph{Bioinformatics}. doi: \href{https://doi.org/10.1093/bioinformatics/btz852}{10.1093/bioinformatics/btz852}.}
#' }
#'
#' @return A \code{\link[ggplot2]{ggplot}} object for the genotype plot.
#'
#' @export
plot.flexdog <- function(x, use_colorblind = TRUE, ...) {
  assertthat::assert_that(is.flexdog(x))
  if (x$input$model == "s1") {
    p1geno <- x$par$pgeno
    p2geno <- NULL
  } else if (x$input$model == "f1") {
    p1geno <- x$par$p1geno
    p2geno <- x$par$p2geno
  } else {
    p1geno <- NULL
    p2geno <- NULL
  }

  pl <- plot_geno(refvec         = x$input$refvec,
                  sizevec        = x$input$sizevec,
                  ploidy         = x$input$ploidy,
                  geno           = x$geno,
                  seq            = x$seq,
                  bias           = x$bias,
                  maxpostprob    = x$maxpostprob,
                  p1ref          = x$input$p1ref,
                  p1size         = x$input$p1size,
                  p2ref          = x$input$p2ref,
                  p2size         = x$input$p2size,
                  p1geno         = p1geno,
                  p2geno         = p2geno,
                  use_colorblind = use_colorblind)
  return(pl)
}

#' Tests if an argument is a \code{flexdog} object.
#'
#' @param x Anything.
#'
#' @return A logical. \code{TRUE} if \code{x} is a \code{flexdog} object, and \code{FALSE} otherwise.
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' is.flexdog("anything")
#' # FALSE
#'
is.flexdog <- function(x) {
  inherits(x, "flexdog")
}


#' Initialize \code{pivec} for \code{\link{flexdog}} EM algorithm.
#'
#' The key idea here is choosing the pi's so that the two modes
#' have equal probability.
#'
#' @inheritParams flexdog_full
#'
#' @seealso \code{\link{flexdog}} for where this is used.
#'
#' @return A vector of numerics. The initial value of \code{pivec}
#'     used in \code{\link{flexdog_full}}.
#'
#' @keywords internal
#'
#' @author David Gerard
initialize_pivec <- function(ploidy,
                             mode,
                             model = c("hw", "bb", "norm", "ash",
                                       "f1", "s1", "f1pp", "s1pp",
                                       "f1ppdr", "s1ppdr",
                                       "flex", "uniform", "custom")) {
  assertthat::are_equal(1, length(ploidy), length(mode))
  assertthat::are_equal(ploidy %% 1, 0)

  model <- match.arg(model)
  if (model == "flex" | model == "uniform") {
    pivec <- rep(x = 1 / (ploidy + 1), length = ploidy + 1)
  } else if (model == "ash") {
    init_type <- "bin"
    if (init_type == "equi") { ## equimodal
      floor_mode <- floor(mode)
      ceil_mode  <- ceiling(mode)
      d <- sum(1 / (ceil_mode:ploidy - ceil_mode + 1)) /
        (sum(1 / (floor_mode - 0:floor_mode + 1)))
      second_half <- 1 / (ploidy - ceil_mode + 1 + (floor_mode + 1) * d)
      first_half <- second_half * d
      ## assertthat::are_equal(first_half, (1 - (ploidy - ceil_mode + 1) * second_half) / (floor_mode + 1))
      pivec <- c(rep(first_half, length = floor_mode + 1),
                 rep(second_half, length = ploidy - ceil_mode + 1))
    } else if (init_type == "bin") { ## binomial
      pvec_init <- stats::dbinom(x = 0:ploidy, size = ploidy,
                                 prob = floor(mode) / ploidy)
      pivec <- get_uni_rep(pvec_init)$pivec + 10 ^-6
      pivec <- pivec / sum(pivec)
    } else if (init_type == "unif") {
      stopifnot((mode > 0) & (mode < ploidy))

      pvec_init <- rep(1 / (ploidy + 1), length = ploidy + 1)
      floor_mode <- max(floor(mode), 0)
      if (floor_mode > ploidy) {
        floor_mode <- ploidy
      }
      pvec_init[floor_mode + 1] <- pvec_init[floor_mode + 1] + 10^-2
      pvec_init <- pvec_init / sum(pvec_init)
      pivec <- get_uni_rep(pvec_init)$pivec
    } else if (init_type == "piunif") {
      pivec <- rep(1 / (ploidy + 1), length = ploidy + 1)
    }
  } else if (model == "hw" | model == "f1" | model == "s1" | model == "bb" | model == "norm" | model == "f1pp" | model == "s1pp" | model == "f1ppdr" | model == "s1ppdr" | model == "custom") {
    if (mode < 0 | mode > 1) {
      stop('initialize_pivec: when model = "hw", mode should be between 0 and 1.\n It is the initialization of the allele frequency.')
    }
    pivec <- stats::dbinom(x = 0:ploidy, size = ploidy, prob = mode)
  } else {
    stop("initialize_pivec: How did you get here?")
  }
  return(pivec)
}


#' Update the distribution of genotypes from various models.
#'
#' @param weight_vec \code{colSums(wik_mat)} from \code{\link{flexdog}}.
#'     This is the sum of current posterior probabilities of each individual
#'     having genotype k.
#' @param model What model are we assuming? See the description in \code{\link{flexdog}} for details.
#' @param control A list of anything else needed to be passed.
#'     E.g. if \code{model = "ash"},
#'     then \code{inner_weights} needs to be passed through \code{control}
#'     (see \code{\link{get_inner_weights}} for how to get this matrix).
#'
#' @return A list with the following elements
#' \describe{
#'   \item{\code{pivec}}{The estimate of the genotype distribution.}
#'   \item{\code{par}}{A list of estimated parameters. An empty list if the model does not contain any parameters other than \code{pivec}.}
#' }
#'
#' @keywords internal
#'
#' @author David Gerard
flex_update_pivec <- function(weight_vec,
                              model = c("hw", "bb", "norm",
                                        "ash", "f1", "s1",
                                        "f1pp", "s1pp",
                                        "f1ppdr", "s1ppdr",
                                        "flex", "uniform", "custom"),
                              control) {
  ## Check input -------------------------------
  ploidy <- length(weight_vec) - 1
  model <- match.arg(model)
  assertthat::assert_that(is.list(control))

  ## Get pivec ---------------------------------
  return_list <- list()
  if (model == "flex") {
    pivec <- weight_vec / sum(weight_vec)
    return_list$pivec <- pivec
    return_list$par <- list()
  } else if (model == "ash") {
    if (is.null(control$inner_weights)) {
      stop('flex_update_pivec: control$inner_weights cannot be NULL when model = "ash"')
    }
    pivec <- uni_em(weight_vec = weight_vec,
                    lmat       = control$inner_weights,
                    pi_init    = control$pivec,
                    itermax    = 200,
                    obj_tol    = 10 ^ -4,
                    lambda     = control$lambda)
    return_list$pivec <- pivec
    return_list$par <- list()
  } else if (model == "hw") {
    alpha <- sum(0:ploidy * weight_vec) / (ploidy * sum(weight_vec))
    pivec <- stats::dbinom(x = 0:ploidy, size = ploidy, prob = alpha)
    return_list$pivec     <- pivec
    return_list$par       <- list()
    return_list$par$alpha <- alpha
  } else if (model == "f1") {
    optim_best       <- list()
    optim_best$value <- -Inf
    for (i in 0:ploidy) { ## parent 1
      for (j in 0:ploidy) { ## parent 2
        pvec <- control$qarray[i + 1, j + 1, ]
        optim_out <- list() ## called this for historical reasons.
        optim_out$value <- f1_obj(alpha = control$fs1_alpha,
                                  pvec = pvec,
                                  weight_vec = weight_vec)
        optim_out$par   <- control$fs1_alpha
        if (!is.null(control$p1_lbb)) {
          optim_out$value <- optim_out$value + control$p1_lbb[i + 1]
        }
        if (!is.null(control$p2_lbb)) {
          optim_out$value <- optim_out$value + control$p2_lbb[j + 1]
        }
        if (optim_out$value > optim_best$value) {
          optim_best <- optim_out
          optim_best$ell1 <- i
          optim_best$ell2 <- j
        }
      }
    }
    return_list$pivec <- (1 - optim_best$par) * control$qarray[optim_best$ell1 + 1, optim_best$ell2 + 1, ] +
      optim_best$par / (ploidy + 1)
    return_list$par <- list()
    return_list$par$p1geno <- optim_best$ell1
    return_list$par$p2geno <- optim_best$ell2
    return_list$par$alpha  <- optim_best$par
    if (control$outliers) {
      return_list$par$out_prop <- control$weight_out / (control$weight_out + sum(weight_vec))
    }
  } else if (model == "s1") {
    optim_best       <- list()
    optim_best$value <- -Inf
    for (i in 0:ploidy) { ## parent
      pvec <- control$qarray[i + 1, i + 1, ]
      optim_out <- list() ## named this for historical reasons
      optim_out$value <- f1_obj(alpha = control$fs1_alpha,
                                pvec = pvec,
                                weight_vec = weight_vec)
      optim_out$par   <- control$fs1_alpha
      if (!is.null(control$p1_lbb)) {
        optim_out$value <- optim_out$value + control$p1_lbb[i + 1]
      }
      if (optim_out$value > optim_best$value) {
        optim_best <- optim_out
        optim_best$ell <- i
      }
    }
    return_list$pivec <- (1 - optim_best$par) * control$qarray[optim_best$ell + 1, optim_best$ell + 1, ] +
      optim_best$par / (ploidy + 1)
    return_list$par <- list()
    return_list$par$pgeno <- optim_best$ell
    return_list$par$alpha  <- optim_best$par
    if (control$outliers) {
      return_list$par$out_prop <- control$weight_out / (control$weight_out + sum(weight_vec))
    }
  } else if (model == "uniform") {
    return_list$pivec <- rep(x = 1 / (ploidy + 1), length = ploidy + 1)
    return_list$par <- list()
  } else if (model == "bb") {
    boundary_val <- 10 ^ -8
    optim_out <- stats::optim(par        = c(control$alpha, control$tau),
                              fn         = obj_for_weighted_lbb,
                              gr         = grad_for_weighted_lbb,
                              method     = "L-BFGS-B",
                              lower      = c(boundary_val, boundary_val),
                              upper      = c(1 - boundary_val, 1 - boundary_val),
                              weight_vec = weight_vec,
                              ploidy     = ploidy,
                              control = list(fnscale = -1))

    return_list$pivec <- dbetabinom(x    = 0:ploidy,
                                    size = ploidy,
                                    mu   = optim_out$par[1],
                                    rho  = optim_out$par[2],
                                    log  = FALSE)
    return_list$par <- list()
    return_list$par$alpha <- optim_out$par[1]
    return_list$par$tau   <- optim_out$par[2]
  } else if (model == "norm") {
    optim_out <- stats::optim(par        = c(control$mu, control$sigma),
                              fn         = obj_for_weighted_lnorm,
                              gr         = grad_for_weighted_lnorm,
                              method     = "L-BFGS-B",
                              lower      = c(-1, 10 ^ -8),
                              upper      = c(ploidy + 1, Inf),
                              weight_vec = weight_vec,
                              ploidy     = ploidy,
                              control    = list(fnscale = -1))
    return_list$pivec <- stats::dnorm(x = 0:ploidy,
                                      mean = optim_out$par[1],
                                      sd = optim_out$par[2],
                                      log = TRUE)
    return_list$pivec <- exp(return_list$pivec - log_sum_exp(return_list$pivec))
    return_list$par       <- list()
    return_list$par$mu    <- optim_out$par[1]
    return_list$par$sigma <- optim_out$par[2]
  } else if (model == "f1pp") {
    temp_list <- update_pp_f1(weight_vec = weight_vec,
                              control    = control)
    return_list$par                 <- list()
    return_list$pivec               <- temp_list$pivec
    return_list$par$p1geno          <- temp_list$p1geno
    return_list$par$p1_pair_weights <- temp_list$p1_pair_weights[[temp_list$p1geno + 1]]
    return_list$par$p2geno          <- temp_list$p2geno
    return_list$par$p2_pair_weights <- temp_list$p2_pair_weights[[temp_list$p2geno + 1]]
  } else if (model == "s1pp") {
    temp_list <- update_pp_s1(weight_vec = weight_vec,
                              control    = control)
    return_list$par                 <- list()
    return_list$pivec               <- temp_list$pivec
    return_list$par$p1geno          <- temp_list$p1geno
    return_list$par$p1_pair_weights <- temp_list$p1_pair_weights[[temp_list$p1geno + 1]]
  } else if (model == "f1ppdr" | model == "s1ppdr") {
    temp_list <- update_dr(weight_vec = weight_vec,
                           model      = model,
                           control    = control)
    return_list$par                  <- list()
    return_list$pivec                <- temp_list$pivec
    return_list$par$p1geno           <- temp_list$p1geno
    return_list$par$p1_pair_weights  <- temp_list$p1_pair_weights[[temp_list$p1geno + 1]]
    return_list$par$p1_dr_proportion <- sum(temp_list$p1_pair_weights[[temp_list$p1geno + 1]][!control$blist$penvec[control$blist$lvec == temp_list$p1geno]])
    if (model == "f1ppdr") {
      return_list$par$p2geno           <- temp_list$p2geno
      return_list$par$p2_pair_weights  <- temp_list$p2_pair_weights[[temp_list$p2geno + 1]]
      return_list$par$p2_dr_proportion <- sum(temp_list$p2_pair_weights[[temp_list$p2geno + 1]][!control$blist$penvec[control$blist$lvec == temp_list$p2geno]])
    }
  } else if (model == "custom") {
    return_list$pivec <- control$prior_vec
    return_list$par <- list()
  } else {
    stop("flex_update_pivec: how did you get here?")
  }
  return(return_list)
}

#' Get the representation of a discrete unimodal probability distribution.
#'
#' NB: In \code{get_uni_rep}, we count the mode starting at 1.
#'     In \code{\link{get_probk_vec}}, we count the mode starting at 0.
#'
#' @param probvec A probability vector. It assumes the probabilities
#'     are ordered according to the ordering of the discrete set.
#'
#' @seealso \code{\link{get_probk_vec}} with option \code{model = "ash"}
#'     for the inverse of this function.
#'
#' @return A list with the following elements
#' \describe{
#' \item{\code{pivec}}{The mixing weights for the unimodal representation.}
#' \item{\code{mode}}{The central value of the unimodal distribution.}
#' }
#'
#' @keywords internal
#'
#' @author David Gerard
#'
get_uni_rep <- function(probvec) {
  assertthat::are_equal(sum(probvec), 1)
  assertthat::assert_that(all(probvec >= 0))

  mode <- which.max(probvec)
  n    <- length(probvec)

  if (mode < n) {
    fak_vec <- c(mode:1, seq_len(n - mode))
  } else {
    fak_vec <- n:1
  }

  if (mode < n) {
    aug_probvec <- c(0, probvec, 0)
    diffvec <- c(aug_probvec[2:(mode + 1)] - aug_probvec[seq_len(mode)],
                 aug_probvec[(mode + 2):(n + 1)] - aug_probvec[(mode + 3):(n + 2)])
  } else {
    aug_probvec <- c(0, probvec)
    diffvec <- aug_probvec[2:(n + 1)] - aug_probvec[seq_len(n)]
  }

  return_list <- list()
  return_list$pivec <- diffvec * fak_vec
  return_list$mode = mode + 0.5
  return(return_list)
}



#' Penalty on pivec used when \code{model = "ash"} in \code{\link{flexdog}}.
#'
#' @param lambda The penalty.
#' @param pivec The vector of mixing proportions for the component
#'     discrete uniform distributions.
#'
#' @return A penalty on the ash mixing weights.
#'
#' @author David Gerard
#'
#' @keywords internal
#'
ashpen_fun <- function(lambda, pivec) {
  lambda * sum(log(pivec))
}

#' Return the probabilities of an offspring's genotype given its
#' parental genotypes for all possible combinations of parental and
#' offspring genotypes. This is for species with polysomal inheritance
#' and bivalent, non-preferential pairing.
#'
#' @param ploidy A positive integer. The ploidy of the species.
#'
#' @author David Gerard
#'
#' @return An three-way array of proportions. The (i, j, k)th element
#'     is the probability of an offspring having k - 1 reference
#'     alleles given that parent 1 has i - 1 reference alleles and
#'     parent 2 has j - 1 reference alleles. Each dimension of the
#'     array is \code{ploidy + 1}. In the dimension names, "A" stands
#'     for the reference allele and "a" stands for the
#'     alternative allele.
#'
#' @examples
#' qarray <- get_q_array(6)
#' apply(qarray, c(1, 2), sum) ## should all be 1's.
#'
#' @export
#'
get_q_array <- function(ploidy) {
  assertthat::assert_that(ploidy > 0)
  assertthat::are_equal(ploidy %% 2, 0)

  qarray <- array(0, dim = rep(ploidy + 1, 3))

  for(oindex in 0:ploidy) {
    for (p1index in 0:ploidy) {
      for (p2index in 0:ploidy) {
        if (p1index + p2index < oindex) {
          qarray[p1index + 1, p2index + 1, oindex + 1] <- 0
        } else {
          minval <- max(0, oindex - p2index)
          maxval <- min(ploidy / 2, p1index)
          aseq <- minval:maxval

          p1prob <- stats::dhyper(x = aseq, m = p1index, n = ploidy - p1index, k = ploidy / 2)
          p2prob <- stats::dhyper(x = oindex - aseq, m = p2index, n = ploidy - p2index, k = ploidy / 2)
          qarray[p1index + 1, p2index + 1, oindex + 1] <- sum(p1prob * p2prob)
        }
      }
    }
  }

  ## get dimnames

  dimvec <- get_dimname(ploidy)
  dimnames(qarray) <- list(parent1 = dimvec, parent2 = dimvec, offspring = dimvec)


  assertthat::assert_that(all(abs(apply(qarray, c(1, 2), sum) - 1) < 10 ^ -14))

  return(qarray)
}


#' Returns a vector character strings that are all of the possible
#' combinations of the reference allele and the non-reference allele.
#'
#' @param ploidy The ploidy of the species.
#'
#' @return For example, if \code{ploidy = 3} then this will return
#'     c("aaa", "Aaa", "AAa", "AAA")
#'
#' @author David Gerard
#'
#' @keywords internal
#'
#'
get_dimname <- function(ploidy) {
  dimvec <- vapply(mapply(FUN = c, lapply(X = 0:ploidy, FUN = rep.int, x = "A"),
                          lapply(X = ploidy:0, FUN = rep.int, x = "a"),
                          SIMPLIFY = FALSE),
                   FUN = paste, collapse = "", FUN.VALUE = "character")
  return(dimvec)
}


