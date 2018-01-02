### plotting for mupdog object




#' Draw a genotype plot from the output of \code{\link{mupdog}}.
#'
#' @param x A mupdog object.
#' @param index The column number of the gene to plot.
#' @param ... Not used.
#'
#' @export
#'
#' @author David Gerard
plot.mupdog <- function(x, index, ...) {
  assertthat::assert_that(is.mupdog(x))

  pl <- updog::plot_geno(ocounts = x$input$refmat[, index], osize = x$input$sizemat[, index],
                         ploidy = x$input$ploidy, ogeno = x$map_dosage[, index],
                         seq_error = x$seq[index], bias_val = x$bias[index],
                         prob_ok = x$maxpostprob[, index], use_colorblind = FALSE) +
    ggplot2::guides(alpha=ggplot2::guide_legend(title="maxpostprob"))
  return(pl)
}


#' Provide some summaries from the output of \code{\link{mupdog}}.
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
