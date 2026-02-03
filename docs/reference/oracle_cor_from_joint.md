# Calculate the correlation of the oracle estimator with the true genotype from the joint distribution matrix.

Calculates the correlation between the oracle MAP estimator (where we
have perfect knowledge about the data generation process) and the true
genotype. This is a useful approximation when you have a lot of
individuals.

## Usage

``` r
oracle_cor_from_joint(jd)
```

## Arguments

- jd:

  A matrix of numerics. Element (i, j) is the probability of genotype
  i - 1 and estimated genotype j - 1. This is usually obtained from
  [`oracle_joint`](https://dcgerard.github.io/updog/reference/oracle_joint.md).

## Value

The Pearson correlation between the true genotype and the oracle
estimator.

## References

- Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018).
  Genotyping Polyploids from Messy Sequencing Data. *Genetics*, 210(3),
  789-807.
  [doi:10.1534/genetics.118.301468](https://doi.org/10.1534/genetics.118.301468)
  .

## See also

[`oracle_joint`](https://dcgerard.github.io/updog/reference/oracle_joint.md)
for getting `jd`.
[`oracle_cor`](https://dcgerard.github.io/updog/reference/oracle_cor.md)
for not having to first calculate `jd`.

## Author

David Gerard

## Examples

``` r
## Hardy-Weinberg population with allele-frequency of 0.75.
## Moderate bias and moderate overdispersion.
ploidy <- 6
dist <- stats::dbinom(0:ploidy, ploidy, 0.75)
jd <- oracle_joint(n = 100, ploidy = ploidy, seq = 0.001,
                   bias = 0.7, od = 0.01, dist = dist)
oracle_cor_from_joint(jd = jd)
#> [1] 0.940216

## Compare to oracle_cor
oracle_cor(n = 100, ploidy = ploidy, seq = 0.001,
           bias = 0.7, od = 0.01, dist = dist)
#> [1] 0.940216

```
