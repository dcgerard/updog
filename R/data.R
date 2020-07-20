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
#' @references
#' \itemize{
#'   \item{Uitdewilligen, J. G., Wolters, A. M. A., Bjorn, B., Borm, T. J., Visser, R. G., & van Eck, H. J. (2013). \href{https://doi.org/10.1371/journal.pone.0062355}{A next-generation sequencing method for genotyping-by-sequencing of highly heterozygous autotetraploid potato.} \emph{PLoS One}, 8(5), e62355.}
#' }
#'
#' @return A list. See the Format Section.
#'
#' @source \url{https://doi.org/10.1371/journal.pone.0062355}
#'
#'
"uitdewilligen"

#' GBS data from Shirasawa et al (2017)
#'
#' Contains counts of reference alleles and total read counts from the GBS data of Shirasawa et al (2017) for
#' the three SNPs used as examples in Gerard et. al. (2018).
#'
#' @format A \code{tibble} with 419 rows and 4 columns:
#' \describe{
#'     \item{id}{The identification label of the individuals.}
#'     \item{snp}{The SNP label.}
#'     \item{counts}{The number of read-counts that support the reference allele.}
#'     \item{size}{The total number of read-counts at a given SNP.}
#' }
#'
#' @return A \code{tibble}. See the Format Section.
#'
#' @source \url{https://doi.org/10.1038/srep44207}
#'
#' @references
#' \itemize{
#'   \item{Shirasawa, Kenta, Masaru Tanaka, Yasuhiro Takahata, Daifu Ma, Qinghe Cao, Qingchang Liu, Hong Zhai, Sang-Soo Kwak, Jae Cheol Jeong, Ung-Han Yoon, Hyeong-Un Lee, Hideki Hirakawa, and Sahiko Isobe "A high-density SNP genetic map consisting of a complete set of homologous groups in autohexaploid sweetpotato (Ipomoea batatas)." \emph{Scientific Reports 7} (2017). DOI: 10.1038/srep44207}
#'   \item{Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018). Genotyping Polyploids from Messy Sequencing Data. \emph{Genetics}, 210(3), 789-807. doi: \href{https://doi.org/10.1534/genetics.118.301468}{10.1534/genetics.118.301468}.}
#' }
#'
"snpdat"

