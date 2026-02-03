# The joint probability of the genotype and the genotype estimate of an oracle estimator.

This returns the joint distribution of the true genotypes and an oracle
estimator given perfect knowledge of the data generating process. This
is a useful approximation when you have a lot of individuals.

## Usage

``` r
oracle_joint(n, ploidy, seq, bias, od, dist)
```

## Arguments

- n:

  The read-depth.

- ploidy:

  The ploidy of the individual.

- seq:

  The sequencing error rate.

- bias:

  The allele-bias.

- od:

  The overdispersion parameter.

- dist:

  The distribution of the alleles.

## Value

A matrix. Element (i, j) is the joint probability of estimating the
genotype to be i+1 when the true genotype is j+1. That is, the estimated
genotype indexes the rows and the true genotype indexes the columns.
This is when using an oracle estimator.

## Details

To come up with `dist`, you need some additional assumptions. For
example, if the population is in Hardy-Weinberg equilibrium and the
allele frequency is `alpha` then you could calculate `dist` using the R
code: `dbinom(x = 0:ploidy, size = ploidy, prob = alpha)`.
Alternatively, if you know the genotypes of the individual's two parents
are, say, `ref_count1` and `ref_count2`, then you could use the
[`get_q_array`](https://dcgerard.github.io/updog/reference/get_q_array.md)
function from the updog package:
`get_q_array(ploidy)[ref_count1 + 1, ref_count2 + 1, ]`.

See the Examples to see how to reconcile the output of `oracle_joint`
with
[`oracle_mis`](https://dcgerard.github.io/updog/reference/oracle_mis.md)
and
[`oracle_mis_vec`](https://dcgerard.github.io/updog/reference/oracle_mis_vec.md).

## References

- Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018).
  Genotyping Polyploids from Messy Sequencing Data. *Genetics*, 210(3),
  789-807.
  [doi:10.1534/genetics.118.301468](https://doi.org/10.1534/genetics.118.301468)
  .

## See also

- [`oracle_plot`](https://dcgerard.github.io/updog/reference/oracle_plot.md):

  For visualizing the joint distribution output from `oracle_joint`.

- [`oracle_mis_from_joint`](https://dcgerard.github.io/updog/reference/oracle_mis_from_joint.md):

  For obtaining the same results as
  [`oracle_mis`](https://dcgerard.github.io/updog/reference/oracle_mis.md)
  directly from the output of `oracle_joint`.

- [`oracle_mis_vec_from_joint`](https://dcgerard.github.io/updog/reference/oracle_mis_vec_from_joint.md):

  For obtaining the same results as
  [`oracle_mis_vec`](https://dcgerard.github.io/updog/reference/oracle_mis_vec.md)
  directly from the output of `oracle_joint`.

- [`oracle_cor_from_joint`](https://dcgerard.github.io/updog/reference/oracle_cor_from_joint.md):

  For obtaining the same results as
  [`oracle_cor`](https://dcgerard.github.io/updog/reference/oracle_cor.md)
  directly from the output of `oracle_joint`.

## Author

David Gerard

## Examples

``` r
## Hardy-Weinberg population with allele-frequency of 0.75.
## Moderate bias and moderate overdispersion.
ploidy <- 4
dist <- stats::dbinom(0:ploidy, ploidy, 0.75)
jd <- oracle_joint(n = 100, ploidy = ploidy, seq = 0.001,
                   bias = 0.7, od = 0.01, dist = dist)
jd
#>              [,1]         [,2]         [,3]         [,4]         [,5]
#> [1,] 3.905665e-03 1.759022e-07 8.379767e-17 1.784566e-29 3.601827e-52
#> [2,] 5.849980e-07 4.379346e-02 2.180335e-03 2.159655e-09 3.235599e-26
#> [3,] 1.897235e-20 3.081362e-03 1.961102e-01 1.099225e-02 6.173803e-14
#> [4,] 1.314974e-34 2.427440e-09 1.264700e-02 4.105964e-01 2.601245e-04
#> [5,] 2.284980e-57 7.543647e-25 9.090651e-12 2.863373e-04 3.161461e-01

## Get same output as oracle_mis this way:
1 - sum(diag(jd))
#> [1] 0.02944818
oracle_mis(n = 100, ploidy = ploidy, seq = 0.001,
           bias = 0.7, od = 0.01, dist = dist)
#> [1] 0.02944818

## Get same output as oracle_mis_vec this way:
1 - diag(sweep(x = jd, MARGIN = 2, STATS = colSums(jd), FUN = "/"))
#> [1] 0.0001497595 0.0657395175 0.0702925658 0.0267344300 0.0008221220
oracle_mis_vec(n = 100, ploidy = ploidy, seq = 0.001,
               bias = 0.7, od = 0.01, dist = dist)
#> [1] 0.0001497595 0.0657395175 0.0702925658 0.0267344300 0.0008221220
```
