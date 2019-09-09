# updog 1.1.1

- Updates documentation to include Gerard and Ferrão (2019) 
  \<[doi:10.1101/751784](https://doi.org/10.1101/751784)\>
  as a reference.
- Minor fixes to documentation.

# updog 1.1.0

- Introduces more flexible priors for more general populations.
- Places a normal prior distribution on the logit of the
  overdispersion parameter. This might change genotype calls from
  previous versions of updog. To reproduce the genotype calls from
  previous versions of updog, simply set `mean_od = 0` and `var_od =
  Inf` in `flexdog()`.
- Adds the `method = "custom"` option to `flexdog()`. This lets users
  choose the genotype distribution if it is completely known a priori.
- Documentation updates.

# updog 1.0.1

- Fixes a bug with option `model = "s1pp"` in `flexdog()`. I was
  originally not constraining the levels of preferential pairing to be
  the same in both segregations of the same parent. This is now
  fixed. But the downside is that `model = "s1pp"` is now only
  supported for `ploidy = 4` or `ploidy = 6`. This is because the
  optimization becomes more difficult for larger ploidy levels.
- I fixed some documentation. Perhaps the biggest error comes from
  this snippet from the original documentation of `flexdog`:

    > The value of `prop_mis` is a very intuitive measure for the
    > quality of the SNP. `prop_mis` is the posterior proportion of
    > individuals mis-genotyped. So if you want only SNPS that
    > accurately genotype, say, 95% of the individuals, you could
    > discard all SNPs with a `prop_mis` under 0.95.

    This now says

    > The value of prop_mis is a very intuitive measure for the
    > quality of the SNP. prop_mis is the posterior proportion of
    > individuals mis-genotyped. So if you want only SNPS that
    > accurately genotype, say, 95% of the individuals, you could
    > discard all SNPs with a prop_mis **over 0.05**.
- I've now exported some C++ functions that I think are useful. You
  can call them in the usual way:
  <http://r-pkgs.had.co.nz/src.html#cpp-import>.


# updog 0.99.0

- This is a complete re-working of the code in `updog`. The old
  version may be found in the [`updogAlpha`
  package](https://github.com/dcgerard/updogAlpha).
- The main function is now `flexdog()`.
- An experimental approach `mupdog()` is now live. We provide no
  guarantees about `mupdog()`'s performance.
- Oracle misclassification error rates may be calculated in
  `oracle_mis()`.
- Genotypes can be simulated using `rgeno()`.
- Next-generation sequencing data can be simulated using `rflexdog()`.
