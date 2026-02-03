# Calculate oracle misclassification error rate.

Given perfect knowledge of the data generating parameters, `oracle_mis`
calculates the misclassification error rate, where the error rate is
taken over both the data generation process and the allele-distribution.
This is an ideal level of the misclassification error rate and any real
method will have a larger rate than this. This is a useful approximation
when you have a lot of individuals.

## Usage

``` r
oracle_mis(n, ploidy, seq, bias, od, dist)
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

A double. The oracle misclassification error rate.

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

## References

- Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018).
  Genotyping Polyploids from Messy Sequencing Data. *Genetics*, 210(3),
  789-807.
  [doi:10.1534/genetics.118.301468](https://doi.org/10.1534/genetics.118.301468)
  .

## Author

David Gerard

## Examples

``` r
## Hardy-Weinberg population with allele-frequency of 0.75.
## Moderate bias and moderate overdispersion.
## See how oracle misclassification error rates change as we
## increase the ploidy.
ploidy <- 2
dist <- stats::dbinom(0:ploidy, ploidy, 0.75)
oracle_mis(n = 100, ploidy = ploidy, seq = 0.001,
           bias = 0.7, od = 0.01, dist = dist)
#> [1] 1.262647e-06

ploidy <- 4
dist <- stats::dbinom(0:ploidy, ploidy, 0.75)
oracle_mis(n = 100, ploidy = ploidy, seq = 0.001,
           bias = 0.7, od = 0.01, dist = dist)
#> [1] 0.02944818

ploidy <- 6
dist <- stats::dbinom(0:ploidy, ploidy, 0.75)
oracle_mis(n = 100, ploidy = ploidy, seq = 0.001,
           bias = 0.7, od = 0.01, dist = dist)
#> [1] 0.1329197
```
