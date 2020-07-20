### plotting for mupdog object


#' Make a genotype plot.
#'
#' The x-axis is the counts of the non-reference allele,
#' and the y-axis is the counts of the reference allele.
#' Transparency is controlled by the \code{maxpostprob}
#' vector. These types of plots are used in Gerard et. al. (2018) and
#' Gerard and Ferrão (2020).
#'
#' If parental genotypes are provided (\code{p1geno} and \code{p2geno}) then
#' they will be colored the same as the offspring. Since they are often hard to see,
#' a small black dot will also indicate their position.
#'
#' @param refvec A vector of non-negative integers. The number of
#'     reference reads observed in the individuals
#' @param sizevec A vector of positive integers. The total number of
#'     reads in the individuals.
#' @param p1ref A vector of non-negative integers. The number of
#'     reference reads observed in parent 1 (if the individuals are
#'     all siblings).
#' @param p1size A vector of positive integers. The total number of
#'     reads in parent 1 (if the individuals are
#'     all siblings).
#' @param p2ref A vector of non-negative integers. The number of
#'     reference reads observed in parent 2 (if the individuals are
#'     all siblings).
#' @param p2size A vector of positive integers. The total number of
#'     reads in parent 2 (if the individuals are
#'     all siblings).
#' @param ploidy A non-negative integer. The ploidy of the species.
#' @param geno The individual genotypes.
#' @param seq The sequencing error rate.
#' @param bias The bias parameter.
#' @param maxpostprob A vector of the posterior probabilities of being at
#'     the modal genotype.
#' @param p1geno Parent 1's genotype.
#' @param p2geno Parent 2's genotype.
#' @param use_colorblind A logical. Should we use a colorblind safe palette (\code{TRUE}),
#'     or not (\code{FALSE})? Only allowed if \code{ploidy <= 6}.
#'
#' @export
#'
#' @references
#' \itemize{
#'   \item{Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018). Genotyping Polyploids from Messy Sequencing Data. \emph{Genetics}, 210(3), 789-807. doi: \href{https://doi.org/10.1534/genetics.118.301468}{10.1534/genetics.118.301468}.}
#'   \item{Gerard, David, and Luís Felipe Ventorim Ferrão. "Priors for genotyping polyploids." Bioinformatics 36, no. 6 (2020): 1795-1800. \href{https://doi.org/10.1093/bioinformatics/btz852}{DOI:10.1093/bioinformatics/btz852}.}
#' }
#'
#' @author David Gerard
#'
#' @return A \code{\link[ggplot2]{ggplot}} object for the genotype plot.
#'
#' @examples
#' data("snpdat")
#' refvec  <- snpdat$counts[snpdat$snp == "SNP1"]
#' sizevec <- snpdat$size[snpdat$snp == "SNP1"]
#' ploidy  <- 6
#' plot_geno(refvec = refvec, sizevec = sizevec, ploidy = ploidy)
#'
plot_geno <- function(refvec,
                      sizevec,
                      ploidy,
                      p1ref          = NULL,
                      p1size         = NULL,
                      p2ref          = NULL,
                      p2size         = NULL,
                      geno           = NULL,
                      seq            = 0,
                      bias           = 1,
                      maxpostprob    = NULL,
                      p1geno         = NULL,
                      p2geno         = NULL,
                      use_colorblind = TRUE) {

  if (ploidy > 6 & use_colorblind) {
    warning("use_colorblind is only supported when ploidy <= 6.\nSwitching to use_colorblind = FALSE.")
    use_colorblind <- FALSE
  }

  assertthat::assert_that(all(refvec >= 0, na.rm = TRUE))
  assertthat::assert_that(all(sizevec >= refvec, na.rm = TRUE))
  assertthat::assert_that(ploidy >= 1)
  assertthat::assert_that(seq >= 0)

  ## get probabilities
  pk <- xi_fun(p = (0:ploidy) / ploidy, eps = seq, h = bias)

  dfdat <- data.frame(A = refvec, a = sizevec - refvec)
  maxcount <- max(max(dfdat$A, na.rm = TRUE), max(dfdat$a, na.rm = TRUE)) + 1
  if (!is.null(geno)) {
    stopifnot(length(geno) == length(refvec))
    dfdat$genotype <- addNA(factor(geno, levels = 0:ploidy))
  }

  if (!is.null(maxpostprob)) {
    stopifnot(all(maxpostprob >= 0, na.rm = TRUE))
    stopifnot(all(maxpostprob <= 1, na.rm = TRUE))
    dfdat$maxpostprob <- maxpostprob
  }

  slopevec <- pk / (1 - pk)
  xend <- pmin(rep(maxcount, ploidy + 1), maxcount / slopevec)
  yend <- pmin(rep(maxcount, ploidy + 1), maxcount * slopevec)
  df_lines <- data.frame(x = rep(0, ploidy + 1), y = rep(0, ploidy + 1),
                         xend = xend, yend = yend)

  ## Plot children
  if (is.null(maxpostprob)) {
    pl <- ggplot2::ggplot(data = dfdat, mapping = ggplot2::aes_string(y = "A", x = "a"))
  } else {
    pl <- ggplot2::ggplot(data = dfdat, mapping = ggplot2::aes_string(y = "A", x = "a", alpha = "maxpostprob"))
  }

  ## add offspring genotypes ------------------------------------------------
  if (!is.null(geno)) {
    pl <- pl + ggplot2::geom_point(ggplot2::aes_string(color = "genotype"))
  } else {
    pl <- pl + ggplot2::geom_point()
  }

  pl <- pl + ggplot2::theme_bw() +
    ggplot2::xlim(0, maxcount) +
    ggplot2::ylim(0, maxcount) +
    ggplot2::ylab("Counts A") +
    ggplot2::xlab("Counts a")  +
    ggplot2::geom_segment(data = df_lines, mapping = ggplot2::aes_string(x = "x", y = "y", xend = "xend", yend = "yend"),
                          lty = 2, alpha = 1/2, color = "black", size = 0.5)

  ## add parents if we have them
  if (!is.null(p1size) & !is.null(p1ref)) {
    stopifnot(all(p1ref >= 0, na.rm = TRUE))
    stopifnot(all(p1size >= p1ref, na.rm = TRUE))
    if (any(p1ref > maxcount | p1size - p1ref > maxcount, na.rm = TRUE)) {
      bad_parent_points <- p1ref > maxcount | p1size - p1ref > maxcount
      bad_parent_points[is.na(bad_parent_points)] <- FALSE
      parent_mult <- (maxcount - 1) / pmax(p1ref[bad_parent_points], p1size[bad_parent_points] - p1ref[bad_parent_points])
      p1ref[bad_parent_points] <- parent_mult * p1ref[bad_parent_points]
      p1size[bad_parent_points]   <- parent_mult * p1size[bad_parent_points]
      pl <- pl + ggplot2::annotate("text", x = p1size[bad_parent_points] - p1ref[bad_parent_points], y = p1ref[bad_parent_points], label = "(scaled)", alpha = 0.3)
    }
    p1dat <- data.frame(A = p1ref, a = p1size - p1ref)
    if (!is.null(p1geno)) {
      p1dat$genotype <- factor(p1geno, levels = 0:ploidy)
      pl <- pl + ggplot2::geom_point(data = p1dat, mapping = ggplot2::aes_string(color = "genotype"),
                                     size = 3, pch = 3, alpha = 1, show.legend = FALSE)
      pl <- pl + ggplot2::geom_point(data = p1dat, size = 1, color = "black", pch = 16, alpha = 1)
    } else {
      pl <- pl + ggplot2::geom_point(data = p1dat, size = 3, color = "black", pch = 3, alpha = 1)
    }
  }
  if (!is.null(p2size) & !is.null(p2ref)) {
    stopifnot(all(p2ref >= 0, na.rm = TRUE))
    stopifnot(all(p2size >= p2ref, na.rm = TRUE))
    if (any(p2ref > maxcount | p2size - p2ref > maxcount, na.rm = TRUE)) {
      bad_parent_points <- p2ref > maxcount | p2size - p2ref > maxcount
      bad_parent_points[is.na(bad_parent_points)] <- FALSE
      parent_mult <- (maxcount - 1) / pmax(p2ref[bad_parent_points], p2size[bad_parent_points] - p2ref[bad_parent_points])
      p2ref[bad_parent_points] <- parent_mult * p2ref[bad_parent_points]
      p2size[bad_parent_points]   <- parent_mult * p2size[bad_parent_points]
      pl <- pl + ggplot2::annotate("text", x = p2size[bad_parent_points] - p2ref[bad_parent_points], y = p2ref[bad_parent_points], label = "(scaled)", alpha = 0.3)
    }
    p2dat <- data.frame(A = p2ref, a = p2size - p2ref)
    if (!is.null(p2geno)) {
      p2dat$genotype <- factor(p2geno, levels = 0:ploidy)
      pl <- pl + ggplot2::geom_point(data = p2dat, mapping = ggplot2::aes_string(color = "genotype"),
                                     size = 3, pch = 4, alpha = 1, show.legend = FALSE)
      pl <- pl + ggplot2::geom_point(data = p2dat, size = 1, color = "black", pch = 16, alpha = 1)
    } else {
      pl <- pl + ggplot2::geom_point(data = p2dat, size = 3, color = "black", pch = 4, alpha = 1)
    }
  }


  ## Set color scale based on use_colorblind --------------------------------------
  if (!is.null(geno) | !is.null(p1geno) | !is.null(p2geno)) {
    if (use_colorblind & requireNamespace("ggthemes", quietly = TRUE)) {
      possible_colors <- paste0(ggthemes::colorblind_pal()(ploidy + 1))
      possible_colors <- possible_colors[length(possible_colors):1]
      pl <- pl + ggplot2::scale_color_manual(values = possible_colors, drop = FALSE, na.translate = TRUE, na.value = "gray50")
    } else if (use_colorblind) {
      pl <- pl + ggplot2::scale_color_hue(drop = FALSE, na.translate = TRUE, na.value = "gray50")
      warning("ggthemes needs to be installed to set use_colorblind = TRUE.")
    } else {
      pl <- pl + ggplot2::scale_color_hue(drop = FALSE, na.translate = TRUE, na.value = "gray50")
    }
  }

  ## Set transparency scale --------------------------------------------------------
  if (!is.null(maxpostprob)) {
    pl <- pl + ggplot2::scale_alpha_continuous(breaks = c(-0.01, 0.25, 0.5, 0.75, 1.01))
  }

  return(pl)
}



