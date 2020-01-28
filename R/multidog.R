################
## Multiple SNP genotyping
################


#' Fit \code{\link{flexdog}} to multiple SNP's.
#'
#' This is a convenience function that will run \code{\link{flexdog}} over many SNP's.
#' Support is provided for parallel computing through the doParallel package.
#' This function has not been extensively tested. Please report any bugs to
#' \url{http://github.com/dcgerard/updog/issues}.
#'
#' You should format your reference counts and total read counts in two
#' separate matrices. The rows should index the markers (SNP's) and the
#' columns should index the individuals. Row names are how we ID the SNP's
#' and column names are how we ID the individuals, and so they are required
#' attributes.
#'
#' See the details of \code{\link{flexdog}} for the possible values of
#' \code{model}.
#'
#' If \code{model = "f1"}, \code{model = "f1pp"}, \code{model = "s1"},
#' or \code{model = "s1pp"} then the user may provide the individual ID
#' for parent(s) via the \code{p1_id} and \code{p2_id} arguments.
#'
#' The output is a list containing two data frames. The first data frame,
#' called \code{snpdf}, contains information on each SNP, such as the allele bias
#' and the sequencing error rate. The second data frame, called \code{inddf},
#' contains information on each individual at each SNP, such as the estimated
#' genotype and the posterior probability of being classified correctly.
#'
#' Using an \code{nc} value greater than \code{1} will allow you to
#' run \code{\link{flexdog}} in parallel. Only set \code{nc} greater than
#' \code{1} if you are sure you have access to the proper number of cores.
#' The upper bound on the value of \code{nc} you should try can be determined
#' by running \code{parallel::detectCores()} in R.
#'
#' SNP's that contain 0 reads (or all missing data) are entirely removed.
#'
#' @inheritParams flexdog
#' @param refmat A matrix of reference read counts. The columns index
#'     the individuals and the rows index the markers (SNP's). This matrix must have
#'     rownames (for the names of the markers) and column names (for the names
#'     of the individuals). These names must match the names in \code{sizemat}.
#' @param sizemat A matrix of total read counts. The columns index
#'     the individuals and the rows index the markers (SNP's). This matrix must have
#'     rownames (for the names of the markers) and column names (for the names
#'     of the individuals). These names must match the names in \code{refmat}.
#' @param nc The number of computing cores to use. This should never be
#'     more than the number of cores available in your computing environment.
#'     You can determine the maximum number of available cores by running
#'     \code{parallel::detectCores()} in R.
#' @param p1_id The ID of the first parent. This should be a character of
#'     length 1. This should correspond to a single column name in \code{refmat}
#'     and \code{sizemat}.
#' @param p2_id The ID of the second parent. This should be a character of
#'     length 1. This should correspond to a single column name in \code{refmat}
#'     and \code{sizemat}.
#'
#' @return A list-like object of two data frames.
#' \describe{
#' \item{\code{snpdf}}{A data frame containing properties of the SNP's (markers).
#'     The rows index the SNP's. The variables include:
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
#'     }}
#' }
#'
#' @author David Gerard
#'
#' @examples
#' \dontrun{
#' data("uitdewilligen")
#' mout <- multidog(refmat = t(uitdewilligen$refmat),
#'                  sizemat = t(uitdewilligen$sizemat),
#'                  ploidy = uitdewilligen$ploidy,
#'                  nc = 2)
#' mout$inddf
#' mout$snpdf
#' }
#'
#' @export
multidog <- function(refmat,
                     sizemat,
                     ploidy,
                     model = c("norm", "hw", "bb", "ash", "s1", "s1pp",
                               "f1", "f1pp", "flex", "uniform", "custom"),
                     nc = 1,
                     p1_id = NULL,
                     p2_id = NULL,
                     bias_init = exp(c(-1, -0.5, 0, 0.5, 1)),
                     outliers = FALSE,
                     prior_vec = NULL,
                     ...) {
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
  assertthat::assert_that(is.numeric(nc))
  nc <- round(nc)
  assertthat::assert_that(nc >= 1)
  assertthat::assert_that(length(nc) == 1)
  if (!is.null(p1_id)) {
    stopifnot(is.character(p1_id))
    stopifnot(length(p1_id) == 1)
    stopifnot(is.element(el = p1_id, set = colnames(refmat)))
    stopifnot(model == "f1" | model == "f1pp" | model == "s1" | model == "s1pp")
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

  ## Get list of SNP's ---------------------------------------------------------
  snplist <- rownames(refmat)

  ## Register workers ----------------------------------------------------------
  if (nc == 1) {
    foreach::registerDoSEQ()
  } else {
    cl = parallel::makeCluster(nc)
    doParallel::registerDoParallel(cl = cl)
    if (foreach::getDoParWorkers() == 1) {
      stop("multidog: nc > 1 but only one core registered from foreach::getDoParWorkers().")
    }
  }

  ## Fit flexdog on all SNP's --------------------------------------------------
  i <- 1
  outlist <- foreach::foreach(i = seq_along(snplist),
                              .export = c("flexdog", "get_bivalent_probs")) %dopar% {
                                current_snp <- snplist[[i]]

                                refvec <- refmat[current_snp, indlist, drop = TRUE]
                                sizevec <- sizemat[current_snp, indlist, drop = TRUE]

                                if (!is.null(p1_id)) {
                                  p1_ref <- refmat[current_snp, p1_id, drop = TRUE]
                                  p1_size <- sizemat[current_snp, p1_id, drop = TRUE]
                                } else {
                                  p1_ref <- NULL
                                  p1_size <- NULL
                                }

                                if (!is.null(p2_id)) {
                                    p2_ref <- refmat[current_snp, p2_id, drop = TRUE]
                                    p2_size <- sizemat[current_snp, p2_id, drop = TRUE]
                                } else {
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
                                                outliers  = outliers,
                                                prior_vec = prior_vec,
                                                ...)


                                names(fout$gene_dist)  <- paste0("Pr_", seq(0, ploidy, by = 1))
                                colnames(fout$postmat) <- paste0("Pr_", seq(0, ploidy, by = 1))

                                snpprop <- cbind(
                                  data.frame(snp      = current_snp,
                                             bias     = fout$bias,
                                             seq      = fout$seq,
                                             od       = fout$od,
                                             prop_mis = fout$prop_mis,
                                             num_iter = fout$num_iter,
                                             llike    = fout$llike,
                                             ploidy   = fout$input$ploidy,
                                             model    = fout$input$model),
                                  as.data.frame(matrix(fout$gene_dist, nrow = 1, dimnames = list(NULL, names(fout$gene_dist))))
                                )


                                if (model == "f1pp") {
                                  blist <- get_bivalent_probs(fout$input$ploidy)

                                  p1weightvec <- rep(0, length = length(blist$lvec))
                                  names(p1weightvec) <- paste0("p1(", apply(blist$pmat, 1, paste, collapse = ","), ")")
                                  p1weightvec[blist$lvec == fout$par$p1geno] <- fout$par$p1_pair_weights

                                  p2weightvec <- rep(0, length = length(blist$lvec))
                                  names(p2weightvec) <- paste0("p2(", apply(blist$pmat, 1, paste, collapse = ","), ")")
                                  p2weightvec[blist$lvec == fout$par$p2geno] <- fout$par$p2_pair_weights

                                  fout$par$p1_pair_weights <- NULL
                                  fout$par$p2_pair_weights <- NULL

                                  par_vec_output <- c(unlist(fout$par), p1weightvec, p2weightvec)
                                  snpprop <- cbind(snpprop, as.data.frame(matrix(par_vec_output, nrow = 1, dimnames = list(NULL, names(par_vec_output)))))

                                } else if (model == "s1pp") {
                                  blist <- get_bivalent_probs(fout$input$ploidy)

                                  p1weightvec <- rep(0, length = length(blist$lvec))
                                  names(p1weightvec) <- paste0("p1(", apply(blist$pmat, 1, paste, collapse = ","), ")")
                                  p1weightvec[blist$lvec == fout$par$p1geno] <- fout$par$p1_pair_weights

                                  fout$par$p1_pair_weights <- NULL

                                  par_vec_output <- c(unlist(fout$par), p1weightvec)
                                  snpprop <- cbind(snpprop, as.data.frame(matrix(par_vec_output, nrow = 1, dimnames = list(NULL, names(par_vec_output)))))

                                } else if (length(fout$par) > 0) {
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
                                  fout$postmat)

                                list(indprop, snpprop)
                              }
  if (nc > 1) {
    parallel::stopCluster(cl)
  }

  inddf <- do.call("rbind", lapply(outlist, function(x) x[[1]]))
  snpdf <- do.call("rbind", lapply(outlist, function(x) x[[2]]))

  retlist <- list(snpdf = as.data.frame(snpdf),
                  inddf = as.data.frame(inddf))
  class(retlist) <- "multidog"
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
#' select which SNP's to plot.
#'
#' @param x The output of \code{\link{multidog}}.
#' @param indices A vector of integers. The indices of the SNP's to plot.
#' @param ... not used.
#'
#' @author David Gerard
#'
#' @export
#'
#' @seealso \code{\link{plot_geno}}.
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


    pllist[[i]] <- plot_geno(refvec         = refvec,
                             sizevec        = sizevec,
                             ploidy         = ploidy,
                             geno           = geno,
                             seq            = seq,
                             bias           = bias,
                             maxpostprob    = maxpostprob,
                             use_colorblind = ploidy <= 6) +
      ggplot2::ggtitle(current_snp)
  }

  return(pllist)
}


#' Return matricized elements from the output of \code{\link{multidog}}.
#'
#' This function will allow you to have genotype estimates, maximum posterior
#' probability, and other values in the form of a matrix.
#'
#' @param x The output of \code{multidog}.
#' @param varname The variable whose values populate the cells. This should
#'     be a single columne from \code{x$inddf}.
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
  assertthat::are_equal(length(varname), 1)
  assertthat::assert_that(varname %in% colnames(x$inddf))

  matout <- reshape2::acast(data      = x$inddf[, c("ind", "snp", varname)],
                            formula   = snp ~ ind,
                            value.var = varname)

  return(matout)
}
