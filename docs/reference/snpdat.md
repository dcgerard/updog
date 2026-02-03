# GBS data from Shirasawa et al (2017)

Contains counts of reference alleles and total read counts from the GBS
data of Shirasawa et al (2017) for the three SNPs used as examples in
Gerard et. al. (2018).

## Usage

``` r
snpdat
```

## Format

A `tibble` with 419 rows and 4 columns:

- id:

  The identification label of the individuals.

- snp:

  The SNP label.

- counts:

  The number of read-counts that support the reference allele.

- size:

  The total number of read-counts at a given SNP.

## Source

[doi:10.1038/srep44207](https://doi.org/10.1038/srep44207)

## Value

A `tibble`. See the Format Section.

## References

- Shirasawa, Kenta, Masaru Tanaka, Yasuhiro Takahata, Daifu Ma, Qinghe
  Cao, Qingchang Liu, Hong Zhai, Sang-Soo Kwak, Jae Cheol Jeong, Ung-Han
  Yoon, Hyeong-Un Lee, Hideki Hirakawa, and Sahiko Isobe "A high-density
  SNP genetic map consisting of a complete set of homologous groups in
  autohexaploid sweetpotato (Ipomoea batatas)." *Scientific Reports 7*
  (2017). [doi:10.1038/srep44207](https://doi.org/10.1038/srep44207)

- Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018).
  Genotyping Polyploids from Messy Sequencing Data. *Genetics*, 210(3),
  789-807.
  [doi:10.1534/genetics.118.301468](https://doi.org/10.1534/genetics.118.301468)
  .
