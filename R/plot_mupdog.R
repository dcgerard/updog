### plotting for mupdog object




#' Draw a genotype plot from the output of \code{\link{mupdog}}.
#'
#' A wrapper for \code{\link[updog]{plot_geno}}. This will create a genotype plot for a single SNP.
#'
#' On a genotype plot, the x-axis contains the counts of the non-reference allele and the y-axis
#' contains the counts of the reference allele. The dashed lines are the expected counts (both reference and alternative)
#' given the sequencing error rate and the allele-bias. The plots are color-coded by the maximum-a-posterior genotypes.
#' Transparency is proportional to the maximum posterior probability for an
#' individual's genotype. Thus, we are less certain of the genotype of more transparent individuals.
#'
#' @param x A mupdog object.
#' @param index The column number of the gene to plot.
#' @param ... Not used.
#'
#' @export
#'
#' @seealso
#' \describe{
#'   \item{\code{\link[updog]{plot_geno}}}{The underlying plotting function.}
#'   \item{\code{\link{mupdog}}}{Creates a \code{mupdog} object.}
#' }
#'
#' @author David Gerard
#'
#' @examples
#' data(mupout)
#' plot(mupout, 4)
#'
plot.mupdog <- function(x, index, ...) {
  assertthat::assert_that(is.mupdog(x))

  pl <- updog::plot_geno(ocounts = x$input$refmat[, index], osize = x$input$sizemat[, index],
                         ploidy = x$input$ploidy, ogeno = x$map_dosage[, index],
                         seq_error = x$seq[index], bias_val = x$bias[index],
                         prob_ok = x$maxpostprob[, index], use_colorblind = FALSE) +
    ggplot2::guides(alpha=ggplot2::guide_legend(title="maxpostprob"))
  return(pl)
}


#' Provides some summaries from the output of \code{\link{mupdog}}.
#'
#' Returns mean-dosage for each individual at each SNP and the SNP-specific distribution of MAP dosages.
#' The mean-dosage in particular might be of interest for downstream analyses.
#'
#' @param object A mupdog object.
#' @param ... Not used.
#'
#' @return A list with the following elements:
#' \itemize{
#' \item \code{freq} A matrix with the frequency of the dosages (rows) for each SNP (columns).
#' \item \code{mean_dosage} The posterior mean dosage of an individual (rows) for each SNP (columns).
#' }
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' data(mupout)
#' msum <- summary(mupout)
#' msum$freq[, 1:5]
#' boxplot(msum$mean_dosage ~ mupout$map_dosage,
#'         xlab = "MAP Dosage", ylab = "Mean Dosage")
#'
summary.mupdog <- function(object, ...) {
  assertthat::assert_that(is.mupdog(object))
  p           <- ncol(object$input$refmat)
  freq_mat    <- apply(rbind(object$map_dosage, matrix(rep(0:object$input$ploidy, times = p), nrow = object$input$ploidy + 1)), 2, table) - 1
  mean_dosage <- apply(sweep(object$postprob, MARGIN = 3, STATS = 0:object$input$ploidy, FUN = "*"), c(1, 2), sum)
  return_list             <- list()
  return_list$freq        <- freq_mat
  return_list$mean_dosage <- mean_dosage
  return(return_list)
}
