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
#' @return A list. See the Format Section.
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
#' @return See the Format Section.
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

#' GBS data from Shirasawa et al (2017)
#'
#' Contains counts of reference alleles and total read counts from the GBS data of Shirasawa et al (2017) for
#' the three SNP's used as examples in Gerard, Ferrao, and Stephens (2017).
#'
#' @format A \code{\link[tibble]{tibble}} with 419 rows and 4 columns:
#' \describe{
#'     \item{id}{The identification label of the individuals.}
#'     \item{snp}{The SNP label.}
#'     \item{counts}{The number of read-counts that support the reference allele.}
#'     \item{size}{The total number of read-counts at a given SNP.}
#' }
#' 
#' @return A \code{\link[tibble]{tibble}}. See the Format Section.
#'
#' @source \url{http://sweetpotato-garden.kazusa.or.jp/}
#'
#' @references Shirasawa, Kenta, Masaru Tanaka, Yasuhiro Takahata, Daifu Ma, Qinghe Cao, Qingchang Liu, Hong Zhai et al. "A high-density SNP genetic map consisting of a complete set of homologous groups in autohexaploid sweetpotato (Ipomoea batatas)." Scientific Reports 7 (2017). DOI: 10.1038/srep44207
#'
#'   Gerard, David, Luis Felipe Ventorim Ferrão, and Matthew Stephens. 2017. "Harnessing Empirical Bayes and Mendelian Segregation for Genotyping Autopolyploids with Messy Sequencing Data." Overleaf Preprint.
#'
"snpdat"

