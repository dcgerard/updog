# updog 2.0.2

This is a massive edit of the updog software. Major changes include:

1. No more support of `model = "ash"`. It seemed that `model = "norm"`
   was always better and faster, so I just got rid of the `"ash"` option.
   This also extremely simplified the code.
2. Removal of `mupdog()`. I think this was a good idea, but the 
   computation was way too slow to be usable.
3. Revision of `model = "f1pp"` and `model = "s1pp"`. These now include
   interpretable parameterizations that are meant to be identified
   via another R package. But support is only for tetraploids right now.
4. `multidog()` now prints some nice ASCII art when it's run.
5. `format_multidog()` now allows you to format multiple variables in terms of a multidimensional array.
6. Fixes a bug where `format_multidog()` was reordering the SNP dimensions. This was fine as long as folks used dimnames properly, but now it should allow folks to also use dim positions.
7. Updog now returns genotype log-likelihoods.


# updog 1.2.1

- Adds `filter_snp()` for filtering the output of `multidog()` based
  on predicates in terms of the variables in `snpdf`.
- Removes stringr from Imports. I was only using it in one place so I
  replaced that code with base R code.
- Removes Rmpfr from Suggests. No longer needed since CVXR is no longer
  suggested.

# updog 1.2.0

- Adds `multidog()` for genotypying multiple SNPs using parallel computing.
- Adds `plot.multidog()` for plotting the output of `multidog()`.
- Adds `format_multidog()` for formatting the output of `multidog()` to be
  a matrix.
- Removes dependency on CVXR. This makes install and maintenance a little easier. 
  The defaults for this specific problem were a little faster anyway.
- No longer changes the color scale in `plot_geno()` based on what 
  genotypes are present.
- In .cpp files, we now coerce objects to be unsigned before comparing. This
  gets rid of some warnings during install.

# updog 1.1.3

- Updates documentation to include the Bioinformatics publication,
  Gerard and Ferrão (2020) 
  \<[doi:10.1093/bioinformatics/btz852](https://doi.org/10.1093/bioinformatics/btz852)\>.
- Adds the "internal" keyword to functions that most users don't need.
- Removes the tidyverse from the Suggests field. I was only using this in 
  the vignettes, so I changed it to base R (except for ggplot2).

# updog 1.1.1

- Updates documentation to include Gerard and Ferrão (2020) 
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
