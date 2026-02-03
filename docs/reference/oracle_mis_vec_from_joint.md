# Get the oracle misclassification error rates (conditional on true genotype) directly from the joint distribution of the genotype and the oracle estimator.

Get the oracle misclassification error rates (conditional on true
genotype) directly from the joint distribution of the genotype and the
oracle estimator.

## Usage

``` r
oracle_mis_vec_from_joint(jd)
```

## Arguments

- jd:

  A matrix of numerics. Element (i, j) is the probability of genotype
  i - 1 and estimated genotype j - 1. This is usually obtained from
  [`oracle_joint`](https://dcgerard.github.io/updog/reference/oracle_joint.md).

## Value

A vector of numerics. Element i is the oracle misclassification error
rate when genotyping an individual with actual genotype i + 1.

## References

- Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018).
  Genotyping Polyploids from Messy Sequencing Data. *Genetics*, 210(3),
  789-807.
  [doi:10.1534/genetics.118.301468](https://doi.org/10.1534/genetics.118.301468)
  .

## See also

[`oracle_joint`](https://dcgerard.github.io/updog/reference/oracle_joint.md)
for getting `jd`.
[`oracle_mis_vec`](https://dcgerard.github.io/updog/reference/oracle_mis_vec.md)
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
oracle_mis_vec_from_joint(jd = jd)
#> [1] 0.001855178 0.186231038 0.262779904 0.249400633 0.177957888 0.103565813
#> [7] 0.005097110

## Compare to oracle_cor
oracle_mis_vec(n = 100, ploidy = ploidy, seq = 0.001,
               bias = 0.7, od = 0.01, dist = dist)
#> [1] 0.001855178 0.186231038 0.262779904 0.249400633 0.177957888 0.103565813
#> [7] 0.005097110
```
