################
## Multiple SNP genotyping
################


#' Fit \code{\link{flexdog}} to multiple SNPs.
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
#' @param nc The number of computing cores to use.
#' @param p1_id The name of the first parent. This should be a character of 
#'     length 1. This should correspond to a single column name in \code{refmat}
#'     and \code{sizemat}.
#' @param p2_id The name of the second parent. This should be a character of 
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
  }
  if (!is.null(p2_id)) {
    stopifnot(is.character(p2_id))
    stopifnot(length(p2_id) == 1)
    stopifnot(is.element(el = p2_id, set = colnames(refmat)))
  }
  if (!is.null(p2_id) & is.null(p1_id)) {
    warning("setting p1_id to be p2_id and setting p2_id to be NULL.")
    p1_id <- p2_id
    p2_id <- NULL
  }
  
  ## Get list of SNP's and individuals -----------------------------------------
  snplist <- rownames(refmat)
  indlist <- colnames(refmat)
  
  if (!is.null(p1_id)) {
    indlist <- indlist[indlist != p1_id]
  }
  if (!is.null(p2_id)) {
    indlist <- indlist[indlist != p2_id]
  }
  
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
                              .export = "flexdog") %dopar% {
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
                                  as.data.frame(matrix(fout$gene_dist, nrow = 1, dimnames = list(NULL, names(fout$gene_dist)))),
                                  as.data.frame(matrix(unlist(fout$par), nrow = 1, dimnames = list(NULL, names(fout$par))))
                                )
                                
                                indprop <- cbind(
                                  data.frame(snp = current_snp,
                                             ind = indlist,
                                             ref = fout$input$refvec,
                                             size = fout$input$sizevec,
                                             geno = fout$geno,
                                             postmean = fout$postmean,
                                             maxpostprob = fout$maxpostprob),
                                  fout$postmat)
                                
                                list(indprop, snpprop)
                              }
  if (nc > 1) {
    parallel::stopCluster(cl)
  }
  
  inddf <- do.call("rbind", lapply(outlist, function(x) x[[1]]))
  snpdf <- do.call("rbind", lapply(outlist, function(x) x[[2]]))
  
  retlist <- list(snpdf = snpdf, 
                  inddf = inddf)
  class(retlist) <- "multidog"
  return(retlist)
}



#' Plot the output of \code{\link{multidog}}.
#' 
#' Produce genotype plots from the output of \code{\link{multidog}}. You may
#' select which SNP's to plot.
#' 
#' @param x The output of \code{\link{multidog}}.
#' @param indices A vector of integers. The indices of the SNP's to plot.
#' @param wait Should we prompt for a new plot (\code{TRUE}) or not 
#'     (\code{FALSE})?
#' @param ... not used.
#' 
#' @author David Gerard
#' 
#' @export
#' 
#' @seealso \code{\link{plot_geno}}.
#' 
plot.multidog <- function(x, indices = seq_len(nrow(x$snpdf)), wait = TRUE, ...) {
  
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
                             use_colorblind = ploidy <= 6) 
    if (wait) {
      print(pllist[[i]])
      readline(prompt = "Press [enter] to continue")
    }
  }
  
  return(pllist)
}

