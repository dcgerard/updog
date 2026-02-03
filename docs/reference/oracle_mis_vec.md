# Returns the oracle misclassification rates for each genotype.

Given perfect knowledge of the data generating parameters,
`oracle_mis_vec` calculates the misclassification error rate at each
genotype. This differs from
[`oracle_mis`](https://dcgerard.github.io/updog/reference/oracle_mis.md)
in that this will *not* average over the genotype distribution to get an
overall misclassification error rate. That is, `oracle_mis_vec` returns
a vector of misclassification error rates *conditional* on each
genotype.

## Usage

``` r
oracle_mis_vec(n, ploidy, seq, bias, od, dist)
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

A vector of numerics. Element i is the oracle misclassification error
rate when genotyping an individual with actual genotype i + 1.

## Details

This is an ideal level of the misclassification error rate and any real
method will have a larger rate than this. This is a useful approximation
when you have a lot of individuals.

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
ploidy <- 4
dist <- stats::dbinom(0:ploidy, ploidy, 0.75)
om <- oracle_mis_vec(n = 100, ploidy = ploidy, seq = 0.001,
                     bias = 0.7, od = 0.01, dist = dist)
om
#> [1] 0.0001497595 0.0657395175 0.0702925658 0.0267344300 0.0008221220

## Get same output as oracle_mis this way:
sum(dist * om)
#> [1] 0.02944818
oracle_mis(n = 100, ploidy = ploidy, seq = 0.001,
           bias = 0.7, od = 0.01, dist = dist)
#> [1] 0.02944818
```
