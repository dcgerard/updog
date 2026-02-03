# Simulate individual genotypes from one of the supported [`flexdog`](https://dcgerard.github.io/updog/reference/flexdog.md) models.

This will simulate genotypes of a sample of individuals drawn from one
of the populations supported by
[`flexdog`](https://dcgerard.github.io/updog/reference/flexdog.md). See
the details of
[`flexdog`](https://dcgerard.github.io/updog/reference/flexdog.md) for
the models allowed. These genotype distributions are described in detail
in Gerard and Ferrão (2020).

## Usage

``` r
rgeno(
  n,
  ploidy,
  model = c("hw", "bb", "norm", "f1", "s1", "flex", "uniform"),
  allele_freq = NULL,
  od = NULL,
  p1geno = NULL,
  p2geno = NULL,
  pivec = NULL,
  mu = NULL,
  sigma = NULL
)
```

## Arguments

- n:

  The number of observations.

- ploidy:

  The ploidy of the species.

- model:

  What form should the prior take? See Details in
  [`flexdog`](https://dcgerard.github.io/updog/reference/flexdog.md).

- allele_freq:

  If `model = "hw"`, then this is the allele frequency of the
  population. For any other model, this should be `NULL`.

- od:

  If `model = "bb"`, then this is the overdispersion parameter of the
  beta-binomial distribution. See
  [`betabinom`](https://dcgerard.github.io/updog/reference/betabinom.md)
  for details. For any other model, this should be `NULL`.

- p1geno:

  Either the first parent's genotype if `model = "f1"`, or the only
  parent's genotype if `model = "s1"`. For any other model, this should
  be `NULL`.

- p2geno:

  The second parent's genotype if `model = "f1"`. For any other model,
  this should be `NULL`.

- pivec:

  A vector of probabilities. If `model = "ash"`, then this represents
  the mixing proportions of the discrete uniforms. If `model = "flex"`,
  then element `i` is the probability of genotype `i - 1`. For any other
  model, this should be `NULL`.

- mu:

  If `model = "norm"`, this is the mean of the normal. For any other
  model, this should be `NULL`.

- sigma:

  If `model = "norm"`, this is the standard deviation of the normal. For
  any other model, this should be `NULL`.

## Value

A vector of length `n` with the genotypes of the sampled individuals.

## Details

List of non-`NULL` arguments:

- `model = "flex"`::

  `pivec`

- `model = "hw"`::

  `allele_freq`

- `model = "f1"`::

  `p1geno` and `p2geno`

- `model = "s1"`::

  `p1geno`

- `model = "uniform"`::

  no non-`NULL` arguments

- `model = "bb"`::

  `allele_freq` and `od`

- `model == "norm"`::

  `mu` and `sigma`

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

## Author

David Gerard

## Examples

``` r
## F1 Population where parent 1 has 1 copy of the referenc allele
## and parent 2 has 4 copies of the reference allele.
ploidy <- 6
rgeno(n = 10, ploidy = ploidy, model = "f1", p1geno = 1, p2geno = 4)
#>  [1] 3 4 1 3 3 3 2 1 4 3

## A population in Hardy-Weinberge equilibrium with an
## allele frequency of 0.75
rgeno(n = 10, ploidy = ploidy, model = "hw", allele_freq = 0.75)
#>  [1] 3 5 5 4 4 4 5 4 6 5
```
