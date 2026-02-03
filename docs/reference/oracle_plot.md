# Construct an oracle plot from the output of [`oracle_joint`](https://dcgerard.github.io/updog/reference/oracle_joint.md).

After obtaining the joint distribution of the true genotype with the
estimated genotype from the oracle estimator using
[`oracle_joint`](https://dcgerard.github.io/updog/reference/oracle_joint.md),
you can use `oracle_plot` to visualize this joint distribution.

## Usage

``` r
oracle_plot(jd)
```

## Arguments

- jd:

  A matrix containing the joint distribution of the true genotype and
  the oracle estimator. Usually, this is obtained by a call from
  [`oracle_joint`](https://dcgerard.github.io/updog/reference/oracle_joint.md).

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html) object
containing the oracle plot. The x-axis indexes the possible values of
the estimated genotype. The y-axis indexes the possible values of the
true genotype. The number in cell (i, j) is the probability that an
individual will have true genotype i but is estimated to have genotype
j. This is when using an oracle estimator. The cells are also
color-coded by the size of the probability in each cell. At the top are
listed the oracle misclassification error rate and the correlation of
the true genotype with the estimated genotype. Both of these quantities
may be derived from the joint distribution.

## References

- Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018).
  Genotyping Polyploids from Messy Sequencing Data. *Genetics*, 210(3),
  789-807.
  [doi:10.1534/genetics.118.301468](https://doi.org/10.1534/genetics.118.301468)
  .

## See also

[`oracle_joint`](https://dcgerard.github.io/updog/reference/oracle_joint.md)
for obtaining `jd`.

## Author

David Gerard

## Examples

``` r
ploidy <- 6
dist <- stats::dbinom(0:ploidy, ploidy, 0.75)
jd <- oracle_joint(n = 100, ploidy = ploidy, seq = 0.001,
                   bias = 0.7, od = 0.01, dist = dist)
pl <- oracle_plot(jd = jd)
print(pl)

```
