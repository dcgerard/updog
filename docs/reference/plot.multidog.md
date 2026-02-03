# Plot the output of [`multidog`](https://dcgerard.github.io/updog/reference/multidog.md).

Produce genotype plots from the output of
[`multidog`](https://dcgerard.github.io/updog/reference/multidog.md).
You may select which SNPs to plot.

## Usage

``` r
# S3 method for class 'multidog'
plot(x, indices = seq(1, min(5, nrow(x$snpdf))), ...)
```

## Arguments

- x:

  The output of
  [`multidog`](https://dcgerard.github.io/updog/reference/multidog.md).

- indices:

  A vector of integers. The indices of the SNPs to plot.

- ...:

  not used.

## Details

On a genotype plot, the x-axis contains the counts of the non-reference
allele and the y-axis contains the counts of the reference allele. The
dashed lines are the expected counts (both reference and alternative)
given the sequencing error rate and the allele-bias. The plots are
color-coded by the maximum-a-posterior genotypes. Transparency is
proportional to the maximum posterior probability for an individual's
genotype. Thus, we are less certain of the genotype of more transparent
individuals. These types of plots are used in Gerard et. al. (2018) and
Gerard and Ferrão (2020).

## References

- Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018).
  Genotyping Polyploids from Messy Sequencing Data. *Genetics*, 210(3),
  789-807.
  [doi:10.1534/genetics.118.301468](https://doi.org/10.1534/genetics.118.301468)
  .

- Gerard, David, and Luís Felipe Ventorim Ferrão. "Priors for genotyping
  polyploids." Bioinformatics 36, no. 6 (2020): 1795-1800.
  [doi:10.1093/bioinformatics/btz852](https://doi.org/10.1093/bioinformatics/btz852)
  .

## See also

[`plot_geno`](https://dcgerard.github.io/updog/reference/plot_geno.md).

## Author

David Gerard
