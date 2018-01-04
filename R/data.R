#' Subset of individuals and SNPs from Uitdewilligen et al (2013).
#'
#' A list containing a matrix of reference counts, a matrix of total counts, and the ploidy level (4) of the species.
#' This is a subset of the data from Uitdewilligen et al (2013).
#'
#' @format A list containing three objects. Two matrices and a numeric scalar:
#' \describe{
#'   \item{refmat}{A matrix of read counts containing the reference allele. The rows index the individuals and the columns index the SNPs.}
#'   \item{sizemat}{A matrix of the total number of read counts. The rows index the individuals and the columns index the SNPs.}
#'   \item{ploidy}{The ploidy level of the species (just 4).}
#' }
#'
#' @references Uitdewilligen, J. G., Wolters, A. M. A., Bjorn, B., Borm, T. J., Visser, R. G., & van Eck, H. J. (2013). \href{https://doi.org/10.1371/journal.pone.0062355}{A next-generation sequencing method for genotyping-by-sequencing of highly heterozygous autotetraploid potato.} PLoS One, 8(5), e62355.
#'
#'
#' @source \url{https://doi.org/10.1371/journal.pone.0062355}
#'
#' @seealso \code{\link{mupout}}: a mupdog fit of these data.
#'
"uitdewilligen"


#' A mupdog fit of the \code{\link{uitdewilligen}} data.
#'
#' @format An object of class \code{\link{mupdog}}.
#'
#' @source The raw data that this was fit to can be found in \code{\link{uitdewilligen}}.
#'
#' @seealso
#' \describe{
#' \item{\code{\link{uitdewilligen}}}{The raw data.}
#' \item{\code{\link{plot.mupdog}}}{A method to plot a \code{\link{mupdog}} object.}
#' \item{\code{\link{summary.mupdog}}}{Calculate some summaries of a \code{\link{mupdog}} object.}
#' \item{\code{\link{mupdog}}}{Function used to create this \code{\link{mupdog}} object.}
#' }
#'
"mupout"
