################
## Multiple SNP genotyping
################


#' Function to combine the parallel output from \code{\link{multidog}()}
#'
#' @param ... Each element should be a list with two data frames. The first
#'     data frames should all have the same column names, and the second
#'     data frames should all have the same column names.
#'
#' @author David Gerard
#'
#' @noRd
combine_flex <- function(...) {
  list(
    do.call(rbind, lapply(list(...), function(x) x[[1]])),
    do.call(rbind, lapply(list(...), function(x) x[[2]]))
  )
}

#' Fit \code{\link{flexdog}} to multiple SNPs.
#'
#' This is a convenience function that will run \code{\link{flexdog}} over many SNPs.
#' Support is provided for parallel computing through the doParallel package.
#' This function has not been extensively tested. Please report any bugs to
#' \url{https://github.com/dcgerard/updog/issues}.
#'
#' You should format your reference counts and total read counts in two
#' separate matrices. The rows should index the markers (SNPs) and the
#' columns should index the individuals. Row names are how we ID the SNPs
#' and column names are how we ID the individuals, and so they are required
#' attributes.
#'
#' If your data are in VCF files, I would recommend importing them using the
#' VariantAnnotation package from Bioconductor
#' \url{https://bioconductor.org/packages/VariantAnnotation/}. It's a great
#' VCF parser.
#'
#' See the details of \code{\link{flexdog}} for the possible values of
#' \code{model}.
#'
#' If \code{model = "f1"}, \code{model = "s1"}, \code{model = "f1pp"}
#' or \code{model = "s1pp"} then the user may
#' provide the individual ID for parent(s) via the \code{p1_id}
#' and \code{p2_id} arguments.
#'
#' The output is a list containing two data frames. The first data frame,
#' called \code{snpdf}, contains information on each SNP, such as the allele bias
#' and the sequencing error rate. The second data frame, called \code{inddf},
#' contains information on each individual at each SNP, such as the estimated
#' genotype and the posterior probability of being classified correctly.
#'
#' SNPs that contain 0 reads (or all missing data) are entirely removed.
#'
#' @section Parallel Computation:
#'
#' The \code{multidog()} function supports parallel computing. It does
#' so through the \href{https://cran.r-project.org/package=future}{future}
#' package.
#'
#' If you are just running \code{multidog()} on a local machine, then you
#' can use the \code{nc} argument to specify the parallelization. Any value
#' of \code{nc} greater than 1 will result in multiple background R sessions to
#' genotype all of the SNPs. The maximum value of \code{nc} you should
#' try can be found by running \code{future::availableCores()}. Running
#' \code{multidog()} using \code{nc} is equivalent to setting the future
#' plan with \code{future::plan(future::multisession, workers = nc)}.
#'
#' Using the future package means that different evaluation strategies
#' are possible. In particular, if you are using a high performance machine,
#' you can explore using the
#' \href{https://cran.r-project.org/package=future.batchtools}{future.batchtools}
#' package to evaluate \code{multidog()} using schedulers like Slurm
#' or TORQUE/PBS.
#'
#' To use a different strategy, set \code{nc = NA} and then
#'  run \code{future::plan()} prior to
#' running \code{multidog()}. For example, to set up forked R processes
#' on your current machine (instead of using background R sessions), you would
#' run (will not work on Windows):
#' \code{future::plan(future::multicore)}, followed by
#' running \code{multidog()} with \code{nc = NA}. See the examples below.
#'
#' @inheritParams flexdog
#' @param refmat A matrix of reference read counts. The columns index
#'     the individuals and the rows index the markers (SNPs). This matrix must have
#'     rownames (for the names of the markers) and column names (for the names
#'     of the individuals). These names must match the names in \code{sizemat}.
#' @param sizemat A matrix of total read counts. The columns index
#'     the individuals and the rows index the markers (SNPs). This matrix must have
#'     rownames (for the names of the markers) and column names (for the names
#'     of the individuals). These names must match the names in \code{refmat}.
#' @param nc The number of computing cores to use when doing parallelization
#'     on your local machine. See the section "Parallel Computation" for how
#'     to implement more complicated evaluation strategies using the
#'     \code{future} package.
#'
#'     When you are specifying other evaluation strategies using the
#'     \code{future} package, you should also set \code{nc = NA}.
#'
#'     The value of \code{nc} should never be
#'     more than the number of cores available in your computing environment.
#'     You can determine the maximum number of available cores by running
#'     \code{future::availableCores()} in R.
#' @param p1_id The ID of the first parent. This should be a character of
#'     length 1. This should correspond to a single column name in \code{refmat}
#'     and \code{sizemat}.
#' @param p2_id The ID of the second parent. This should be a character of
#'     length 1. This should correspond to a single column name in \code{refmat}
#'     and \code{sizemat}.
#'
#' @return A list-like object of two data frames.
#' \describe{
#' \item{\code{snpdf}}{A data frame containing properties of the SNPs (markers).
#'     The rows index the SNPs. The variables include:
#'     \describe{
#'     \item{\code{snp}}{The name of the SNP (marker).}
#'     \item{\code{bias}}{The estimated allele bias of the SNP.}
#'     \item{\code{seq}}{The estimated sequencing error rate of the SNP.}
#'     \item{\code{od}}{The estimated overdispersion parameter of the SNP.}
#'     \item{\code{prop_mis}}{The estimated proportion of individuals
#'         misclassified in the SNP.}
#'     \item{\code{num_iter}}{The number of iterations performed during
#'         the EM algorithm for that SNP.}
#'     \item{\code{llike}}{The maximum marginal likelihood of the SNP.}
#'     \item{\code{ploidy}}{The provided ploidy of the species.}
#'     \item{\code{model}}{The provided model for the prior genotype
#'         distribution.}
#'     \item{\code{p1ref}}{The user-provided reference read counts of parent 1.}
#'     \item{\code{p1size}}{The user-provided total read counts of parent 1.}
#'     \item{\code{p2ref}}{The user-provided reference read counts of parent 2.}
#'     \item{\code{p2size}}{The user-provided total read counts of parent 2.}
#'     \item{\code{Pr_k}}{The estimated frequency of individuals with genotype
#'         k, where k can be any integer between 0 and the ploidy level.}
#'     \item{Model specific parameter estimates}{See the return value of
#'         \code{par} in the help page of \code{\link{flexdog}}.}
#'     }}
#' \item{\code{inddf}}{A data frame containing the properties of the
#'     individuals at each SNP. The variables include:
#'     \describe{
#'     \item{\code{snp}}{The name of the SNP (marker).}
#'     \item{\code{ind}}{The name of the individual.}
#'     \item{\code{ref}}{The provided reference counts for that individual at
#'          that SNP.}
#'     \item{\code{size}}{The provided total counts for that individual at
#'          that SNP.}
#'     \item{\code{geno}}{The posterior mode genotype for that individual at
#'          that SNP. This is the estimated reference allele dosage for a
#'          given individual at a given SNP.}
#'     \item{\code{postmean}}{The posterior mean genotype for that individual
#'          at that SNP. This is a continuous genotype estimate of the
#'          reference allele dosage for a given individual at a given SNP.}
#'     \item{\code{maxpostprob}}{The maximum posterior probability. This
#'          is the posterior probability that the individual was genotyped
#'          correctly.}
#'     \item{\code{Pr_k}}{The posterior probability that a given individual
#'          at a given SNP has genotype k, where k can vary from 0 to the
#'          ploidy level of the species.}
#'     \item{\code{logL_k}}{The genotype \emph{log}-likelihoods for dosage
#'          k for a given individual at a given SNP, where k can vary f
#'          rom 0 to the ploidy level of the species.}
#'     }}
#' }
#'
#' @author David Gerard
#'
#' @seealso
#' \itemize{
#'   \item{\code{\link{flexdog}()}:}{For the underlying genotyping function.}
#'   \item{\code{\link{format_multidog}()}:}{For converting the output
#'       of \code{multidog()} to a matrix.}
#'   \item{\code{\link{filter_snp}()}:}{For filtering SNPs using the
#'       output of \code{multidog()}.}
#' }
#'
#' @examples
#' \dontrun{
#' data("uitdewilligen")
#'
#' ## Run multiple R sessions using the `nc` variable.
#' mout <- multidog(refmat = t(uitdewilligen$refmat),
#'                  sizemat = t(uitdewilligen$sizemat),
#'                  ploidy = uitdewilligen$ploidy,
#'                  nc = 2)
#' mout$inddf
#' mout$snpdf
#'
#' ## Run multiple external R sessions on the local machine.
#' ## Note that we set `nc = NA`.
#' cl <- parallel::makeCluster(2, timeout = 60)
#' future::plan(future::cluster, workers = cl)
#' mout <- multidog(refmat = t(uitdewilligen$refmat),
#'                  sizemat = t(uitdewilligen$sizemat),
#'                  ploidy = uitdewilligen$ploidy,
#'                  nc = NA)
#' mout$inddf
#' mout$snpdf
#'
#' ## Close cluster and reset future to current R process
#' parallel::stopCluster(cl)
#' future::plan(future::sequential)
#' }
#'
#' @export
multidog <- function(refmat,
                     sizemat,
                     ploidy,
                     model = c("norm",
                               "hw",
                               "bb",
                               "s1",
                               "s1pp",
                               "f1",
                               "f1pp",
                               "flex",
                               "uniform",
                               "custom"),
                     nc = 1,
                     p1_id = NULL,
                     p2_id = NULL,
                     bias_init = exp(c(-1, -0.5, 0, 0.5, 1)),
                     prior_vec = NULL,
                     ...) {

  cat(paste0(  "    |                                   *.#,%    ",
             "\n   |||                                 *******/  ",
             "\n |||||||    (**..#**.                  */   **/  ",
             "\n|||||||||    */****************************/*%   ",
             "\n   |||    &****..,*.************************/    ",
             "\n   |||     (....,,,*,...****%********/(******    ",
             "\n   |||                ,,****%////,,,,./.****/    ",
             "\n   |||                  /**//         .*///....  ",
             "\n   |||                  .*/*/%#         .,/   ., ",
             "\n   |||               , **/   #%         .*    .. ",
             "\n   |||                               ,,,*        ",
             "\n\nWorking on it..."))

  ## Check input --------------------------------------------------------------
  assertthat::assert_that(is.matrix(refmat))
  assertthat::assert_that(is.matrix(sizemat))
  assertthat::assert_that(is.numeric(refmat))
  assertthat::assert_that(is.numeric(sizemat))
  assertthat::assert_that(!is.null(colnames(refmat)))
  assertthat::assert_that(!is.null(colnames(sizemat)))
  assertthat::assert_that(!is.null(rownames(refmat)))
  assertthat::assert_that(!is.null(rownames(sizemat)))
  if (length(setdiff(colnames(refmat), colnames(sizemat))) != 0) {
    stop(paste0("multidog: refmat and sizemat must have the same column names\n",
                "The following column names are in the set difference:\n",
                setdiff(colnames(refmat), colnames(sizemat))))
  }
  if (length(setdiff(rownames(refmat), rownames(sizemat))) != 0) {
    stop(paste0("multidog: refmat and sizemat must have the same row names\n",
                "The following row names are in the set difference:\n",
                setdiff(rownames(refmat), rownames(sizemat))))
  }
  model <- match.arg(model)
  assertthat::assert_that(length(nc) == 1)
  if (!is.na(nc)) {
    assertthat::assert_that(is.numeric(nc))
    nc <- round(nc)
    assertthat::assert_that(nc >= 1)
  }
  if (!is.null(p1_id)) {
    stopifnot(is.character(p1_id))
    stopifnot(length(p1_id) == 1)
    stopifnot(is.element(el = p1_id, set = colnames(refmat)))
    stopifnot(model == "f1" | model == "s1" | model == "f1pp" | model == "s1pp")
  }
  if (!is.null(p2_id)) {
    stopifnot(is.character(p2_id))
    stopifnot(length(p2_id) == 1)
    stopifnot(is.element(el = p2_id, set = colnames(refmat)))
    stopifnot(model == "f1" | model == "f1pp")
  }
  if (!is.null(p2_id) & is.null(p1_id)) {
    warning("setting p1_id to be p2_id and setting p2_id to be NULL.")
    p1_id <- p2_id
    p2_id <- NULL
  }

  ## Get list of individuals ---------------------------------------------------
  indlist <- colnames(refmat)

  if (!is.null(p1_id)) {
    indlist <- indlist[indlist != p1_id]
  }
  if (!is.null(p2_id)) {
    indlist <- indlist[indlist != p2_id]
  }

  ## Remove NA SNPs ------------------------------------------------------------
  which_bad_size <- apply(X = (sizemat[, indlist, drop = FALSE] == 0) | is.na(sizemat[, indlist, drop = FALSE]),
                          MARGIN = 1,
                          FUN = all)
  which_bad_ref <- apply(X = is.na(refmat[, indlist, drop = FALSE]),
                         MARGIN = 1,
                         FUN = all)

  bad_snps <- unique(c(rownames(sizemat)[which_bad_size], rownames(refmat)[which_bad_ref]))

  if (length(bad_snps) > 0) {
    if (length(bad_snps) == nrow(sizemat)) {
      stop("multidog: All SNPs are missing.")
    }
    sizemat <- sizemat[!(rownames(sizemat) %in% bad_snps), , drop = FALSE]
    refmat  <- refmat[!(rownames(refmat) %in% bad_snps), , drop = FALSE]
  }

  ## Get list of SNPs ---------------------------------------------------------
  snplist <- rownames(refmat)

  ## Extract parent vectors ---------------------------------------------------
  if (!is.null(p1_id)) {
    p1_refvec <- refmat[, p1_id]
    p1_sizevec <- sizemat[, p1_id]
  } else {
    p1_refvec <- rep(NA_real_, length.out = length(snplist))
    p1_sizevec <- rep(NA_real_, length.out = length(snplist))
  }

  if (!is.null(p2_id)) {
    p2_refvec <- refmat[, p2_id]
    p2_sizevec <- sizemat[, p2_id]
  } else {
    p2_refvec <- rep(NA_real_, length.out = length(snplist))
    p2_sizevec <- rep(NA_real_, length.out = length(snplist))
  }

  refmat <- refmat[, indlist, drop = FALSE]
  sizemat <- sizemat[, indlist, drop = FALSE]

  ## Register doFuture  -------------------------------------------------------
  oldDoPar <- doFuture::registerDoFuture()
  on.exit(with(oldDoPar, foreach::setDoPar(fun=fun, data=data, info=info)), add = TRUE)

  ## Register workers ----------------------------------------------------------
  if (!is.na(nc)) {
    if (nc > 1) {
      oplan <- future::plan(future::multisession, workers = nc)
      on.exit(future::plan(oplan), add = TRUE)
    }
  }

  ## Fit flexdog on all SNPs --------------------------------------------------
  current_snp <- NULL
  refvec <- NULL
  sizevec <- NULL
  p1_ref <- NULL
  p1_size <- NULL
  p2_ref <- NULL
  p2_size <- NULL
  retlist <- foreach::foreach(current_snp   = iterators::iter(snplist),
                              refvec        = iterators::iter(refmat, by = "row"),
                              sizevec       = iterators::iter(sizemat, by = "row"),
                              p1_ref        = iterators::iter(p1_refvec),
                              p1_size       = iterators::iter(p1_sizevec),
                              p2_ref        = iterators::iter(p2_refvec),
                              p2_size       = iterators::iter(p2_sizevec),
                              .export       = c("flexdog"),
                              .combine      = combine_flex,
                              .multicombine = TRUE) %dorng% {

                                if (is.na(p1_ref) || is.na(p1_size)) {
                                  p1_ref <- NULL
                                  p1_size <- NULL
                                }

                                if (is.na(p2_ref) || is.na(p2_size)) {
                                  p2_ref <- NULL
                                  p2_size <- NULL
                                }

                                fout <- flexdog(refvec    = refvec,
                                                sizevec   = sizevec,
                                                ploidy    = ploidy,
                                                model     = model,
                                                p1ref     = p1_ref,
                                                p1size    = p1_size,
                                                p2ref     = p2_ref,
                                                p2size    = p2_size,
                                                snpname   = current_snp,
                                                bias_init = bias_init,
                                                verbose   = FALSE,
                                                prior_vec = prior_vec,
                                                ...
                                                )

                                names(fout$gene_dist)  <- paste0("Pr_", seq(0, ploidy, by = 1))
                                colnames(fout$postmat) <- paste0("Pr_", seq(0, ploidy, by = 1))
                                colnames(fout$genologlike) <- paste0("logL_", seq(0, ploidy, by = 1))

                                ## change to NA so can return in data frame ----
                                if (is.null(p1_ref)) {
                                  p1_ref <- NA_real_
                                }
                                if (is.null(p1_size)) {
                                  p1_size <- NA_real_
                                }
                                if (is.null(p2_ref)) {
                                  p2_ref <- NA_real_
                                }
                                if (is.null(p2_size)) {
                                  p2_size <- NA_real_
                                }

                                snpprop <- cbind(
                                  data.frame(snp      = current_snp,
                                             bias     = fout$bias,
                                             seq      = fout$seq,
                                             od       = fout$od,
                                             prop_mis = fout$prop_mis,
                                             num_iter = fout$num_iter,
                                             llike    = fout$llike,
                                             ploidy   = fout$input$ploidy,
                                             model    = fout$input$model,
                                             p1ref    = p1_ref,
                                             p1size   = p1_size,
                                             p2ref    = p2_ref,
                                             p2size   = p2_size),
                                  as.data.frame(matrix(fout$gene_dist, nrow = 1, dimnames = list(NULL, names(fout$gene_dist))))
                                )


                                if (length(fout$par) > 0) {
                                  par_vec_output <- unlist(fout$par)
                                  snpprop <- cbind(snpprop, as.data.frame(matrix(par_vec_output, nrow = 1, dimnames = list(NULL, names(par_vec_output)))))
                                }


                                indprop <- cbind(
                                  data.frame(snp         = current_snp,
                                             ind         = indlist,
                                             ref         = fout$input$refvec,
                                             size        = fout$input$sizevec,
                                             geno        = fout$geno,
                                             postmean    = fout$postmean,
                                             maxpostprob = fout$maxpostprob),
                                  fout$postmat,
                                  fout$genologlike)

                                list(snpdf = snpprop, inddf = indprop)
                              }

  names(retlist) <- c("snpdf", "inddf")
  attr(retlist, "rng") <- NULL
  attr(retlist, "doRNG_version") <- NULL
  class(retlist) <- "multidog"

  cat("done!")

  return(retlist)
}

#' Tests if an argument is a \code{multidog} object.
#'
#' @param x Anything.
#'
#' @return A logical. \code{TRUE} if \code{x} is a \code{multidog} object, and \code{FALSE} otherwise.
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' is.multidog("anything")
#' # FALSE
#'
is.multidog <- function(x) {
  inherits(x, "multidog")
}



#' Plot the output of \code{\link{multidog}}.
#'
#' Produce genotype plots from the output of \code{\link{multidog}}. You may
#' select which SNPs to plot.
#'
#' On a genotype plot, the x-axis contains the counts of the non-reference allele and the y-axis
#' contains the counts of the reference allele. The dashed lines are the expected counts (both reference and alternative)
#' given the sequencing error rate and the allele-bias. The plots are color-coded by the maximum-a-posterior genotypes.
#' Transparency is proportional to the maximum posterior probability for an
#' individual's genotype. Thus, we are less certain of the genotype of more transparent individuals. These
#' types of plots are used in Gerard et. al. (2018) and Gerard and Ferrão (2020).
#'
#' @param x The output of \code{\link{multidog}}.
#' @param indices A vector of integers. The indices of the SNPs to plot.
#' @param ... not used.
#'
#' @author David Gerard
#'
#' @export
#'
#' @seealso \code{\link{plot_geno}}.
#'
#' @references
#' \itemize{
#'   \item{Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018). Genotyping Polyploids from Messy Sequencing Data. \emph{Genetics}, 210(3), 789-807. \doi{10.1534/genetics.118.301468}.}
#'   \item{Gerard, David, and Luís Felipe Ventorim Ferrão. "Priors for genotyping polyploids." Bioinformatics 36, no. 6 (2020): 1795-1800. \doi{10.1093/bioinformatics/btz852}.}
#' }
#'
plot.multidog <- function(x, indices = seq(1, min(5, nrow(x$snpdf))), ...) {
  assertthat::assert_that(is.multidog(x))

  all(indices <= nrow(x$snpdf))

  pllist <- list()
  for (i in seq_along(indices)) {
    current_index <- indices[[i]]
    current_snp <- x$snpdf$snp[[current_index]]
    bias <- x$snpdf$bias[[current_index]]
    seq <- x$snpdf$seq[[current_index]]
    ploidy <- x$snpdf$ploidy[[current_index]]

    refvec <- x$inddf[x$inddf$snp == current_snp, "ref", drop = TRUE]
    sizevec <- x$inddf[x$inddf$snp == current_snp, "size", drop = TRUE]
    geno <- x$inddf[x$inddf$snp == current_snp, "geno", drop = TRUE]
    maxpostprob <- x$inddf[x$inddf$snp == current_snp, "maxpostprob", drop = TRUE]

    if (is.na(x$snpdf$p1ref[[current_index]]) | is.na(x$snpdf$p1size[[current_index]])) {
      p1_ref <- NULL
      p1_size <- NULL
    } else {
      p1_ref <- x$snpdf$p1ref[[current_index]]
      p1_size <- x$snpdf$p1size[[current_index]]
    }

    if (is.na(x$snpdf$p2ref[[current_index]]) | is.na(x$snpdf$p2size[[current_index]])) {
      p2_ref <- NULL
      p2_size <- NULL
    } else {
      p2_ref <- x$snpdf$p2ref[[current_index]]
      p2_size <- x$snpdf$p2size[[current_index]]
    }

    pllist[[i]] <- plot_geno(refvec         = refvec,
                             sizevec        = sizevec,
                             ploidy         = ploidy,
                             geno           = geno,
                             seq            = seq,
                             bias           = bias,
                             maxpostprob    = maxpostprob,
                             use_colorblind = ploidy <= 6,
                             p1ref          = p1_ref,
                             p1size         = p1_size,
                             p2ref          = p2_ref,
                             p2size         = p2_size) +
      ggplot2::ggtitle(current_snp)
  }

  return(pllist)
}


#' Return arrayicized elements from the output of \code{\link{multidog}}.
#'
#' This function will allow you to have genotype estimates, maximum posterior
#' probability, and other values in the form of a matrix/array. If multiple
#' variable names are provided, the data are formatted as a 3-dimensional
#' array with the dimensions corresponding to (individuals, SNPs, variables).
#'
#' Note that the order of the individuals will be reshuffled. The order of the
#' SNPs should be the same as in \code{x$snpdf}.
#'
#' @param x The output of \code{multidog}.
#' @param varname A character vector of the variable names whose values
#'     populate the cells. These should be column names from \code{x$inddf}.
#'
#' @author David Gerard
#'
#' @export
#'
format_multidog <- function(x, varname = "geno") {
  assertthat::assert_that(is.multidog(x))
  assertthat::assert_that(!is.null(x$inddf))
  assertthat::assert_that(is.data.frame(x$inddf))
  assertthat::assert_that(is.character(varname))
  stopifnot(varname %in% colnames(x$inddf))

  snporder <- x$snpdf[["snp"]]
  if (length(varname) == 1) {
    matout <- reshape2::acast(data      = x$inddf[, c("ind", "snp", varname)],
                              formula   = snp ~ ind,
                              value.var = varname)
    matout <- matout[match(snporder, rownames(matout)), ]
  } else {
    matout <- reshape2::acast(
      data = reshape2::melt(data = x$inddf[, c("ind", "snp", varname)],
                            id.vars = c("ind", "snp"),
                            measure.vars = varname),
      formula = snp ~ ind ~ variable,
      value.var = "value"
    )
    names(dimnames(matout)) <- c("snp", "ind", "variable")
    matout <- matout[match(snporder, dimnames(matout)[["snp"]]), , ]
  }

  return(matout)
}

#' Filter SNPs based on the output of \code{\link{multidog}()}.
#'
#' Filter based on provided logical predicates in terms of the variable
#' names in \code{x$snpdf}. This function filters both \code{x$snpdf}
#' and \code{x$inddf}.
#'
#' @param x The output of \code{multidog}.
#' @param expr Logical predicate expression defined in terms of the variables
#'     in \code{x$snpdf}. Only SNPs where the condition evaluates to
#'     \code{TRUE} are kept.
#'
#' @examples
#' \dontrun{
#' data("uitdewilligen")
#' mout <- multidog(refmat = t(uitdewilligen$refmat),
#'                  sizemat = t(uitdewilligen$sizemat),
#'                  ploidy = uitdewilligen$ploidy,
#'                  nc = 2)
#'
#' ## The following filters are for educational purposes only and should
#' ## not be taken as a default filter:
#' mout2 <- filter_snp(mout, bias < 0.8 & od < 0.003)
#' }
#'
#' @seealso
#' \itemize{
#'   \item{\code{\link{multidog}()}:}{For the variables in \code{x$snpdf}
#'       which you can filter by.}
#' }
#'
#' @author David Gerard
#'
#' @export
filter_snp <- function(x, expr) {
  assertthat::assert_that(is.multidog(x))
  cond <- eval(expr = substitute(expr), envir = x$snpdf)
  x$snpdf <- x$snpdf[cond, , drop = FALSE]
  goodsnps <- x$snpdf$snp
  x$inddf <- x$inddf[x$inddf$snp %in% goodsnps, , drop = FALSE]
  return(x)
}




#' Save the output of \code{\link{multidog}()} to a VCF file.
#'
#' This is an experimental function to save the output of
#' \code{\link{multidog}()} to a VCF file.
#'
#' This function uses the Bioconductor packages VariantAnnotation,
#' GenomicRanges, S4Vectors, and IRanges. You can install these using the
#' BiocManager package via:
#'
#' \code{install.packages("BiocManager")}
#'
#' \code{BiocManager::install(c("VariantAnnotation", "GenomicRanges", "S4Vectors", "IRanges"))}
#'
#' To read more about the VCF format, see the official documentation
#' on the Samtools website: \url{https://samtools.github.io/hts-specs/}.
#'
#' @param obj An object of class \code{\link{multidog}()}.
#' @param filename A string, the path to save the VCF file.
#'
#' @examples
#' \dontrun{
#' library(VariantAnnotation)
#' data("snpdat")
#' refmat <- reshape2::acast(data = snpdat,
#'                           formula = snp ~ id,
#'                           value.var = "counts")
#' sizemat <- reshape2::acast(data = snpdat,
#'                            formula = snp ~ id,
#'                            value.var = "size")
#' mout <- multidog(refmat = refmat,
#'                  sizemat = sizemat,
#'                  ploidy = 6,
#'                  model = "s1",
#'                  p1_id = "Xushu18")
#' export_vcf(obj = mout, filename = "./sweet_potato.vcf")
#' spvcf <- readVcf("./sweet_potato.vcf")
#' }
#'
#'
export_vcf <- function(obj, filename) {
  if (requireNamespace("VariantAnnotation", quietly = TRUE) &&
      requireNamespace("GenomicRanges", quietly = TRUE) &&
      requireNamespace("S4Vectors", quietly = TRUE) &&
      requireNamespace("IRanges", quietly = TRUE)) {

    ploidy <- unique(obj$snpdf$ploidy)
    stopifnot(length(ploidy) == 1)
    obj$inddf$alt <- obj$inddf$size - obj$inddf$ref
    geno <- S4Vectors::SimpleList(
      AD = format_multidog(x = obj, varname = c("alt", "ref")),
      DP = format_multidog(x = obj, varname = "size"),
      DS = format_multidog(x = obj, varname = "postmean"),
      GP = format_multidog(x = obj, varname = paste0("Pr_", 0:ploidy)),
      GL = format_multidog(x = obj, varname = paste0("logL_", 0:ploidy)) / log(10)
    )
    for (i in seq_along(geno)) {
      dimnames(geno[[i]]) <- NULL
    }

    nind <- ncol(geno$AD)
    nsnp <- nrow(obj$snpdf)

    infodf <- S4Vectors::DataFrame(
      row.names = c("snp",
                    "bias",
                    "seq",
                    "od",
                    "prop_mis",
                    "llike",
                    "ploidy",
                    "model",
                    paste0("Pr_", 0:ploidy)),
      Description = c("SNP name",
                      "Allele bias",
                      "Sequencing error rate",
                      "Overdispersion",
                      "Proportion of individuals genotyped incorrectly",
                      "Maximized marginal log-likelihood",
                      "Ploidy",
                      "Model used",
                      paste0("Prior probability of dosage ", 0:ploidy)),
      Type = c("String", "Float", "Float", "Float",
               "Float", "Float", "Integer", "String",
               rep("Float", ploidy + 1)),
      Number = 1)

    vcfobj <- VariantAnnotation::VCF(
      rowRanges = GenomicRanges::GRanges(
        seqnames=paste0("snp:", seq_len(nsnp)),
        ranges=NULL,
        strand=NULL,
        seqinfo = NULL,
        names = obj$snpdf$snp),
      colData = as(matrix(nrow = nind, ncol = 0), "DataFrame"),
      exptData = list(
        header = VariantAnnotation::VCFHeader(
          reference = character(),
          samples = character(),
          header = IRanges::DataFrameList(
            fileformat = S4Vectors::DataFrame(row.names = "fileformat", Value = "VCFv4.3"),
            fileDate = S4Vectors::DataFrame(row.names = "fileDate", Value = gsub("-", "", Sys.Date())),
            source = S4Vectors::DataFrame(row.names = "source", Value = paste0("updogv", utils::packageVersion("updog"))),
            FORMAT = S4Vectors::DataFrame(row.names = c("AD", "DP", "DS", "GP", "GL"),
                                          Number = c("R", "1", "1", "G", "G"),
                                          Type = c("Integer", "Integer", "Float", "Float", "Float"),
                                          Description = c("Read depth for each allele",
                                                          "Read depth",
                                                          "Posterior mean genotype",
                                                          "Genotype posterior probabilities",
                                                          "Genotype likelihoods")),
            INFO = infodf
            )
          )
        ),
      fixed = as(matrix(nrow = nsnp, ncol = 0), "DataFrame"),
      geno = geno,
      info = as(obj$snpdf[, row.names(infodf)], "DataFrame"),
      collapsed = TRUE,
      verbose = TRUE
      )

    VariantAnnotation::writeVcf(obj = vcfobj, filename = filename)
  } else {
    message(paste0("Need to have VariantAnnotation, S4Vectors,\n",
                   "GenomicRanges, and IRanges installed to run",
                   "\nexport_vcf()\n",
                   "To install, run in R:\n\n",
                   "install.packages(\"BiocManager\")\n",
                   "BiocManager::install(c(\"VariantAnnotation\",",
                   " \"GenomicRanges\", \"S4Vectors\", \"IRanges\"))\n"))
  }
}

