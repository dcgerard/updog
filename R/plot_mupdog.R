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
