# Draw a genotype plot from the output of [`flexdog`](https://dcgerard.github.io/updog/reference/flexdog.md).

A wrapper for
[`plot_geno`](https://dcgerard.github.io/updog/reference/plot_geno.md).
This will create a genotype plot for a single SNP.

## Usage

``` r
# S3 method for class 'flexdog'
plot(x, use_colorblind = TRUE, ...)
```

## Arguments

- x:

  A `flexdog` object.

- use_colorblind:

  Should we use a colorblind-safe palette (`TRUE`) or not (`FALSE`)?
  `TRUE` is only allowed if the ploidy is less than or equal to 6.

- ...:

  Not used.

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html) object
for the genotype plot.

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

- [`plot_geno`](https://dcgerard.github.io/updog/reference/plot_geno.md):

  The underlying plotting function.

- [`flexdog`](https://dcgerard.github.io/updog/reference/flexdog.md):

  Creates a `flexdog` object.

## Author

David Gerard
