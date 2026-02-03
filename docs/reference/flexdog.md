# Flexible genotyping for polyploids from next-generation sequencing data.

Genotype polyploid individuals from next generation sequencing (NGS)
data while assuming the genotype distribution is one of several forms.
`flexdog` does this while accounting for allele bias, overdispersion,
sequencing error. The method is described in detail in Gerard et. al.
(2018) and Gerard and Ferrão (2020). See
[`multidog()`](https://dcgerard.github.io/updog/reference/multidog.md)
for running flexdog on multiple SNPs in parallel.

## Usage

``` r
flexdog(
  refvec,
  sizevec,
  ploidy,
  model = c("norm", "hw", "bb", "s1", "s1pp", "f1", "f1pp", "flex", "uniform", "custom"),
  p1ref = NULL,
  p1size = NULL,
  p2ref = NULL,
  p2size = NULL,
  snpname = NULL,
  bias_init = exp(c(-1, -0.5, 0, 0.5, 1)),
  verbose = TRUE,
  prior_vec = NULL,
  ...
)
```

## Arguments

- refvec:

  A vector of counts of reads of the reference allele.

- sizevec:

  A vector of total counts.

- ploidy:

  The ploidy of the species. Assumed to be the same for each individual.

- model:

  What form should the prior (genotype distribution) take? See Details
  for possible values.

- p1ref:

  The reference counts for the first parent if `model = "f1"` or
  `model = "f1pp"`, or for the only parent if `model = "s1"` or
  `model = "s1pp"`.

- p1size:

  The total counts for the first parent if `model = "f1"` or
  `model = "f1pp"`, or for the only parent if `model = "s1"` or
  `model = "s1pp"`.

- p2ref:

  The reference counts for the second parent if `model = "f1"` or
  `model = "f1pp"`.

- p2size:

  The total counts for the second parent if `model = "f1"` or
  `model = "f1pp"`.

- snpname:

  A string. The name of the SNP under consideration. This is just
  returned in the `input` list for your reference.

- bias_init:

  A vector of initial values for the bias parameter over the multiple
  runs of
  [`flexdog_full()`](https://dcgerard.github.io/updog/reference/flexdog_full.md).

- verbose:

  Should we output more (`TRUE`) or less (`FALSE`)?

- prior_vec:

  The pre-specified genotype distribution. Only used if
  `model = "custom"` and must otherwise be `NULL`. If specified, then it
  should be a vector of length `ploidy + 1` with non-negative elements
  that sum to 1.

- ...:

  Additional parameters to pass to
  [`flexdog_full()`](https://dcgerard.github.io/updog/reference/flexdog_full.md).

## Value

An object of class `flexdog`, which consists of a list with some or all
of the following elements:

- `bias`:

  The estimated bias parameter.

- `seq`:

  The estimated sequencing error rate.

- `od`:

  The estimated overdispersion parameter.

- `num_iter`:

  The number of EM iterations ran. You should be wary if this equals
  `itermax`.

- `llike`:

  The maximum marginal log-likelihood.

- `postmat`:

  A matrix of posterior probabilities of each genotype for each
  individual. The rows index the individuals and the columns index the
  allele dosage.

- `genologlike`:

  A matrix of genotype *log*-likelihoods of each genotype for each
  individual. The rows index the individuals and the columns index the
  allele dosage.

- `gene_dist`:

  The estimated genotype distribution. The `i`th element is the
  proportion of individuals with genotype `i-1`.

- `par`:

  A list of the final estimates of the parameters of the genotype
  distribution. The elements included in `par` depends on the value of
  `model`:

  `model = "norm"`:

  :   

      `mu`:

      :   The normal mean.

      `sigma`:

      :   The normal standard deviation (not variance).

  `model = "hw"`:

  :   

      `alpha`:

      :   The major allele frequency.

  `model = "bb"`:

  :   

      `alpha`:

      :   The major allele frequency.

      `tau`:

      :   The overdispersion parameter. See the description of `rho` in
          the Details of
          [`betabinom()`](https://dcgerard.github.io/updog/reference/betabinom.md).

  `model = "s1"`:

  :   

      `pgeno`:

      :   The allele dosage of the parent.

      `alpha`:

      :   The mixture proportion of the discrete uniform (included and
          fixed at a small value mostly for numerical stability
          reasons). See the description of `fs1_alpha` in
          [`flexdog_full()`](https://dcgerard.github.io/updog/reference/flexdog_full.md).

  `model = "f1"`:

  :   

      `p1geno`:

      :   The allele dosage of the first parent.

      `p2geno`:

      :   The allele dosage of the second parent.

      `alpha`:

      :   The mixture proportion of the discrete uniform (included and
          fixed at a small value mostly for numerical stability
          reasons). See the description of `fs1_alpha` in
          [`flexdog_full()`](https://dcgerard.github.io/updog/reference/flexdog_full.md).

  `model = "s1pp"`:

  :   

      `ell1`:

      :   The estimated dosage of the parent.

      `tau1`:

      :   The estimated double reduction parameter of the parent.
          Available if `ell1` is `1`, `2`, or `3`. Identified if `ell1`
          is `1` or `3`.

      `gamma1`:

      :   The estimated preferential pairing parameter. Available if
          `ell1` is `2`. However, it is not returned in an identified
          form.

      `alpha`:

      :   The mixture proportion of the discrete uniform (included and
          fixed at a small value mostly for numerical stability
          reasons). See the description of `fs1_alpha` in
          [`flexdog_full()`](https://dcgerard.github.io/updog/reference/flexdog_full.md).

  `model = "f1pp"`:

  :   

      `ell1`:

      :   The estimated dosage of parent 1.

      `ell2`:

      :   The estimated dosage of parent 2.

      `tau1`:

      :   The estimated double reduction parameter of parent 1.
          Available if `ell1` is `1`, `2`, or `3`. Identified if `ell1`
          is `1` or `3`.

      `tau2`:

      :   The estimated double reduction parameter of parent 2.
          Available if `ell2` is `1`, `2`, or `3`. Identified if `ell2`
          is `1` or `3`.

      `gamma1`:

      :   The estimated preferential pairing parameter of parent 1.
          Available if `ell1` is `2`. However, it is not returned in an
          identified form.

      `gamma2`:

      :   The estimated preferential pairing parameter of parent 2.
          Available if `ell2` is `2`. However, it is not returned in an
          identified form.

      `alpha`:

      :   The mixture proportion of the discrete uniform (included and
          fixed at a small value mostly for numerical stability
          reasons). See the description of `fs1_alpha` in
          [`flexdog_full()`](https://dcgerard.github.io/updog/reference/flexdog_full.md).

  `model = "flex"`:

  :   `par` is an empty list.

  `model = "uniform"`:

  :   `par` is an empty list.

  `model = "custom"`:

  :   `par` is an empty list.

- `geno`:

  The posterior mode genotype. These are your genotype estimates.

- `maxpostprob`:

  The maximum posterior probability. This is equivalent to the posterior
  probability of correctly genotyping each individual.

- `postmean`:

  The posterior mean genotype. In downstream association studies, you
  might want to consider using these estimates.

- `input$refvec`:

  The value of `refvec` provided by the user.

- `input$sizevec`:

  The value of `sizevec` provided by the user.

- `input$ploidy`:

  The value of `ploidy` provided by the user.

- `input$model`:

  The value of `model` provided by the user.

- `input$p1ref`:

  The value of `p1ref` provided by the user.

- `input$p1size`:

  The value of `p1size` provided by the user.

- `input$p2ref`:

  The value of `p2ref` provided by the user.

- `input$p2size`:

  The value of `p2size` provided by the user.

- `input$snpname`:

  The value of `snpname` provided by the user.

- `prop_mis`:

  The posterior proportion of individuals genotyped incorrectly.

## Details

Possible values of the genotype distribution (values of `model`) are:

- `"norm"`:

  A distribution whose genotype frequencies are proportional to the
  density value of a normal with some mean and some standard deviation.
  Unlike the `"bb"` and `"hw"` options, this will allow for
  distributions both more and less dispersed than a binomial. This seems
  to be the most robust to violations in modeling assumptions, and so is
  the default. This prior class was developed in Gerard and Ferrão
  (2020).

- `"hw"`:

  A binomial distribution that results from assuming that the population
  is in Hardy-Weinberg equilibrium (HWE). This actually does pretty well
  even when there are minor to moderate deviations from HWE. Though it
  does not perform as well as the \`"norm"\` option when there are
  severe deviations from HWE.

- `"bb"`:

  A beta-binomial distribution. This is an overdispersed version of
  `"hw"` and can be derived from a special case of the Balding-Nichols
  model.

- `"s1"`:

  This prior assumes the individuals are all full-siblings resulting
  from one generation of selfing. I.e. there is only one parent. This
  model assumes a particular type of meiotic behavior: polysomic
  inheritance with bivalent, non-preferential pairing.

- `"f1"`:

  This prior assumes the individuals are all full-siblings resulting
  from one generation of a bi-parental cross. This model assumes a
  particular type of meiotic behavior: polysomic inheritance with
  bivalent, non-preferential pairing.

- `"f1pp"`:

  This prior allows for double reduction and preferential pairing in an
  F1 population of tretraploids.

- `"s1pp"`:

  This prior allows for double reduction and preferential pairing in an
  S1 population of tretraploids.

- `"flex"`:

  Generically any categorical distribution. Theoretically, this works
  well if you have a lot of individuals. In practice, it seems to be
  much less robust to violations in modeling assumptions.

- `"uniform"`:

  A discrete uniform distribution. This should never be used in
  practice.

- `"custom"`:

  A pre-specified prior distribution. You specify it using the
  `prior_vec` argument. You should almost never use this option in
  practice.

You might think a good default is `model = "uniform"` because it is
somehow an "uninformative prior." But it is very informative and tends
to work horribly in practice. The intuition is that it will estimate the
allele bias and sequencing error rates so that the estimated genotypes
are approximately uniform (since we are assuming that they are
approximately uniform). This will usually result in unintuitive
genotyping since most populations don't have a uniform genotype
distribution. I include it as an option only for completeness. Please
don't use it.

The value of `prop_mis` is a very intuitive measure for the quality of
the SNP. `prop_mis` is the posterior proportion of individuals
mis-genotyped. So if you want only SNPS that accurately genotype, say,
95% of the individuals, you could discard all SNPs with a `prop_mis`
over `0.05`.

The value of `maxpostprob` is a very intuitive measure for the quality
of the genotype estimate of an individual. This is the posterior
probability of correctly genotyping the individual when using `geno`
(the posterior mode) as the genotype estimate. So if you want to
correctly genotype, say, 95% of individuals, you could discard all
individuals with a `maxpostprob` of under `0.95`. However, if you are
just going to impute missing genotypes later, you might consider not
discarding any individuals as `flexdog`'s genotype estimates will
probably be more accurate than other more naive approaches, such as
imputing using the grand mean.

In most datasets I've examined, allelic bias is a major issue. However,
you may fit the model assuming no allelic bias by setting
`update_bias = FALSE` and `bias_init = 1`.

Prior to using `flexdog`, during the read-mapping step, you could try to
get rid of allelic bias by using WASP
([doi:10.1101/011221](https://doi.org/10.1101/011221) ). If you are
successful in removing the allelic bias (because its only source was the
read-mapping step), then setting `update_bias = FALSE` and
`bias_init = 1` would be reasonable. You can visually inspect SNPs for
bias by using
[`plot_geno()`](https://dcgerard.github.io/updog/reference/plot_geno.md).

`flexdog()`, like most methods, is invariant to which allele you label
as the "reference" and which you label as the "alternative". That is, if
you set `refvec` with the number of alternative read-counts, then the
resulting genotype estimates will be the estimated allele dosage of the
alternative allele.

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

Run `browseVignettes(package = "updog")` in R for example usage. Other
useful functions include:

- [`multidog()`](https://dcgerard.github.io/updog/reference/multidog.md):

  For running `flexdog()` on multiple SNPs in parallel.

- [`flexdog_full()`](https://dcgerard.github.io/updog/reference/flexdog_full.md):

  For additional parameter options when running `flexdog()`.

- [`rgeno()`](https://dcgerard.github.io/updog/reference/rgeno.md):

  For simulating genotypes under the allowable prior models in
  `flexdog()`.

- [`rflexdog()`](https://dcgerard.github.io/updog/reference/rflexdog.md):

  For simulating read-counts under the assumed likelihood model in
  `flexdog()`.

- [`plot.flexdog()`](https://dcgerard.github.io/updog/reference/plot.flexdog.md):

  For plotting the output of `flexdog()`.

- [`oracle_mis()`](https://dcgerard.github.io/updog/reference/oracle_mis.md):

  For calculating the oracle genotyping error rates. This is useful for
  read-depth calculations *before* collecting data. After you have data,
  using the value of `prop_mis` is better.

- [`oracle_cor()`](https://dcgerard.github.io/updog/reference/oracle_cor.md):

  For calculating the correlation between the true genotypes and an
  oracle estimator (useful for read-depth calculations *before*
  collecting data).

## Author

David Gerard

## Examples

``` r
# \donttest{

## An S1 population where the first individual
## is the parent.
data("snpdat")
ploidy  <- 6
refvec  <- snpdat$counts[snpdat$snp == "SNP2"]
sizevec <- snpdat$size[snpdat$snp == "SNP2"]
fout    <- flexdog(refvec   = refvec[-1],
                   sizevec  = sizevec[-1],
                   ploidy   = ploidy,
                   model    = "s1",
                   p1ref    = refvec[1],
                   p1size   = sizevec[1])
#>          Fit: 1 of 5 
#> Initial Bias: 0.3678794 
#> Log-Likelihood: -558.9084 
#> Keeping new fit.
#> 
#>          Fit: 2 of 5 
#> Initial Bias: 0.6065307 
#> Log-Likelihood: -519.8545 
#> Keeping new fit.
#> 
#>          Fit: 3 of 5 
#> Initial Bias: 1 
#> Log-Likelihood: -519.8545 
#> Keeping old fit.
#> 
#>          Fit: 4 of 5 
#> Initial Bias: 1.648721 
#> Log-Likelihood: -519.8545 
#> Keeping new fit.
#> 
#>          Fit: 5 of 5 
#> Initial Bias: 2.718282 
#> Log-Likelihood: -519.8545 
#> Keeping new fit.
#> 
#> Done!
plot(fout)
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_point()`).


# }

## A natural population. We will assume a
## normal prior since there are so few
## individuals.
data("uitdewilligen")
ploidy  <- 4
refvec  <- uitdewilligen$refmat[, 1]
sizevec <- uitdewilligen$sizemat[, 1]
fout    <- flexdog(refvec  = refvec,
                   sizevec = sizevec,
                   ploidy  = ploidy,
                   model   = "norm")
#>          Fit: 1 of 5 
#> Initial Bias: 0.3678794 
#> Log-Likelihood: -15.31235 
#> Keeping new fit.
#> 
#>          Fit: 2 of 5 
#> Initial Bias: 0.6065307 
#> Log-Likelihood: -15.31235 
#> Keeping new fit.
#> 
#>          Fit: 3 of 5 
#> Initial Bias: 1 
#> Log-Likelihood: -15.83366 
#> Keeping old fit.
#> 
#>          Fit: 4 of 5 
#> Initial Bias: 1.648721 
#> Log-Likelihood: -15.83359 
#> Keeping old fit.
#> 
#>          Fit: 5 of 5 
#> Initial Bias: 2.718282 
#> Log-Likelihood: -15.83359 
#> Keeping old fit.
#> 
#> Done!
plot(fout)



```
