#' Reference allele count log-likelihood
#'
#' @param gl The genotype log-likelihoods. The rows index the individuals
#'     and the columns index the genotypes.
#'
#' @return A vector of log-likelihoods for the non-reference allele counts at a
#'     single site (so-called SAF likelihood).
#'
#' @examples
#' data("uitdewilligen", package = "updog")
#' mout <- multidog(
#'   refmat = uitdewilligen$refmat,
#'   sizemat = uitdewilligen$sizemat,
#'   ploidy = 4)
#' glarray <- format_multidog(x = mout, varname = paste0("logL_", 0:4))
#' gl <- glarray[1, , ]
#' sfs_ll(gl)
#'
#' @author David Gerard
#'
#' @noRd
sfs_ll <- function(gl) {
  ploidy <- ncol(gl) - 1
  nind <- nrow(gl)
  zmat <- matrix(-Inf, nrow = nind + 1, ncol = nind * ploidy + 1)
  zmat[1, 1] <- 0

  for (j in 1:nind) {
    for (ell in 0:(j * ploidy)) {
      poss <- 0:min(ell, ploidy)
      zmat[j + 1, ell + 1] <- log_sum_exp(zmat[j, ell - poss + 1] + lchoose(n = ploidy, k = poss) + gl[j, poss + 1])
    }
  }

  ll <- zmat[nind + 1, ] - lchoose(n = nind * ploidy, k = 0:(nind * ploidy))
  return(ll)
}

#' Naive calculation of SAF likelihoods
#'
#' Uses Equation (13) of Li (2011). Only use this for samples of
#' size less than 5. This is just used for debugging.
#'
#' @author David Gerard
#'
#' @noRd
saf_naive <- function(gl) {
  ploidy <- ncol(gl) - 1
  nind <- nrow(gl)

  ilist <- rep(list(0:ploidy), nind)
  names(ilist) <- paste0("Ind", 1:nind)
  genodf <- expand.grid(ilist)
  genodf$k <- rowSums(genodf)
  genodf$sumand <- 0
  for (i in seq_len(nrow(genodf))) {
    for (j in seq_len(nind)) {
      gi <- genodf[i, j]
      genodf$sumand[[i]] <- genodf$sumand[[i]] + gl[j, gi + 1] + lchoose(ploidy, gi)
    }
  }

  ll <- rep(NA_real_, nind * ploidy + 1)
  for (k in 0:(nind * ploidy)) {
    ll[[k + 1]] <- log_sum_exp(genodf$sumand[genodf$k == k]) - lchoose(nind * ploidy, k)
  }
  return(ll)
}

#' Reference allele count log-likelihood for all loci
#'
#' @param glarray An array of genotype log-likelihoods. Rows index SNPs,
#'     columns index individuals, faces index genotypes
#'
#' @return A matrix of allele count likelihoods (the so called SAF likelihoods).
#'     The columns index the SNPs and the rows index the allele counts.
#'
#' @examples
#' data("uitdewilligen", package = "updog")
#' mout <- multidog(
#'   refmat = uitdewilligen$refmat,
#'   sizemat = uitdewilligen$sizemat,
#'   ploidy = 4)
#' glarray <- format_multidog(x = mout, varname = paste0("logL_", 0:4))
#' gl <- glarray[1, , ]
#' sfs_ll(gl)
#'
#' @author David Gerard
#'
#' @noRd
sfs_ll_multi <- function(glarray) {
  stopifnot(length(dim(glarray)) == 3)
  nsnps <- dim(glarray)[[1]]
  nind <- dim(glarray)[[2]]
  ploidy <- dim(glarray)[[3]] - 1

  llmat <- matrix(NA_real_, ncol = nsnps, nrow = nind * ploidy + 1)

  for (i in seq_len(nsnps)) {
    llmat[, i] <- sfs_ll(gl = glarray[i, , ])
  }

  return(llmat)
}

#' EM algorithm of Li (2011) Section 2.3.8 to get estimate sfs
#'
#' @param llmat Allele-count log-likelihood matrix from
#'    \code{\link{sfs_ll_multi}()}
#' @param itermax Maximum number of iterations to run em
#' @param tol The stopping tolerance
#' @param verbose Should we print more or less?
#'
#' @author David Gerard
#'
#' @noRd
sfs_em <- function(llmat, itermax = 1000, tol = 1e-5, verbose = FALSE) {
  nsnps <- ncol(llmat)
  phi_log <- rep(log(1 / nrow(llmat)), nrow(llmat)) ## initial SFS

  iter <- 1
  err <- Inf
  while (iter <= itermax && err >= tol) {
    phi_log_old <- phi_log
    phi_log <- apply(apply(phi_log + llmat, 2, function(x) x - log_sum_exp(x)), 1, log_sum_exp) - log(nsnps)

    iter <- iter + 1
    err <- sum(abs(exp(phi_log_old) - exp(phi_log)))

    if (verbose) {
      cat("Iter:", iter - 1, "\n",
          "Err:", err, "\n\n")
    }
  }

  return(exp(phi_log))
}

#' Estimate SFS using algorithm of Li (2011)
#'
#' We use the dynamic algorithm of Section 2.3.5 to obtain the site allele
#' frequency likelihood for all SNPs, then use the EM algorithm of section
#' 2.3.8 to get the site frequency spectrum.
#'
#' @param glarray An array of genotype log-likelihoods. Rows index SNPs,
#'     columns index individuals, faces index genotypes
#' @param tol The stopping tolerance for the EM algorithm
#' @param itermax The maximum number of EM iterations.
#'
#' @return The estimated SFS.
#'
#' @references
#' \itemize{
#'   \item Li, H. (2011). A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics, 27(21), 2987-2993.
#' }
#'
#' @examples
#' \dontrun{
#' data("uitdewilligen", package = "updog")
#' mout <- multidog(
#'   refmat = uitdewilligen$refmat,
#'   sizemat = uitdewilligen$sizemat,
#'   ploidy = 4)
#'
#' ## estimate SFS using genotype likelihoods
#' glarray <- format_multidog(x = mout, varname = paste0("logL_", 0:4))
#' sfs_gl <- sfs_est(glarray)
#'
#' ## compare to SFS using posterior modes
#' genomat <- format_multidog(x = mout, varname = "geno")
#' sfs_mode <- unclass(
#'   prop.table(
#'     table(
#'       factor(
#'         rowSums(genomat),
#'         levels = 0:(ncol(genomat) * 4)
#'       )
#'      )
#'    )
#'  )
#'
#' ## pretty similar, some differences:
#' df <- data.frame(mode = sfs_mode, gl = sfs_gl)
#' round(df[df$mode > 1e-6 | df$gl > 1e-6, ], digits = 3)
#' }
#'
#' @author David Gerard
#'
#' @export
sfs_est <- function(glarray, tol = 1e-5, itermax = 1000) {
  llmat <- sfs_ll_multi(glarray = glarray)
  phihat <- sfs_em(llmat = llmat, tol = tol, itermax = itermax)
  names(phihat) <- 0:(length(phihat) - 1)
  return(phihat)
}
