
<!-- README.md is generated from README.Rmd. Please edit that file -->

# updog

[![Build
Status](https://travis-ci.org/dcgerard/updog.svg?branch=master)](https://travis-ci.org/dcgerard/updog)
[![Appveyor Build
status](https://ci.appveyor.com/api/projects/status/c80fdy61ead6s3vr?svg=true)](https://ci.appveyor.com/project/dcgerard/updog-06s8t)
[![codecov](https://codecov.io/gh/dcgerard/updog/branch/master/graph/badge.svg)](https://codecov.io/gh/dcgerard/updog)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/updog)](https://cran.r-project.org/package=updog)
[![](http://cranlogs.r-pkg.org/badges/grand-total/updog)](https://cran.r-project.org/package=updog)

Updog provides a suite of methods for genotyping polyploids from
next-generation sequencing (NGS) data. It does this while accounting for
many common features of NGS data: allele bias, overdispersion,
sequencing error, and (possibly) outlying observations. It is named
updog for “Using Parental Data for Offspring Genotyping” because we
originally developed the method for full-sib populations, but it works
now for more general populations. The method is described in detail
Gerard et. al. (2018)
\<[doi:10.1534/genetics.118.301468](https://doi.org/10.1534/genetics.118.301468)\>.
Additional details concerning prior specification are described in
Gerard and Ferrão (2019)
\<[doi:10.1093/bioinformatics/btz852](https://doi.org/10.1093/bioinformatics/btz852)\>.

The main function is `flexdog()`, which provides many options for the
distribution of the genotypes in your sample. Novel genotype
distributions include the class of proportional normal distributions
(`model = "norm"`) and the class of discrete unimodal distributions
(`model = "ash"`). The default is `model = "norm"` because it is the
most robust to varying genotype distributions, but feel free to use more
specialized priors if you have more information on the data.

`multidog()` is a convenience function that let’s you run `flexdog()`
over many SNP’s. It has support for parallel computing.

Also provided are:

  - An experimental function `mupdog()`, which allows for correlation
    between the individuals’ genotypes while jointly estimating the
    genotypes of the individuals at all provided SNPs. The
    implementation uses a variational approximation. This is designed
    for samples where the individuals share a complex relatedness
    structure (e.g. siblings, cousins, uncles, half-siblings, etc).
    Right now there are no guarantees about this function’s performance.
  - Functions to simulate genotypes (`rgeno()`) and read-counts
    (`rflexdog()`). These support all of the models available in
    `flexdog()`.
  - Functions to evaluate oracle genotyping performance:
    `oracle_joint()`, `oracle_mis()`, `oracle_mis_vec()`, and
    `oracle_cor()`. We mean “oracle” in the sense that we assume that
    the entire data generation process is known (i.e. the genotype
    distribution, sequencing error rate, allele bias, and overdispersion
    are all known). These are good approximations when there are a lot
    of individuals (but not necessarily large read-depth).

The original `updog` package is now named `updogAlpha` and may be found
[here](https://github.com/dcgerard/updogAlpha).

See also [ebg](https://github.com/pblischak/polyploid-genotyping),
[fitPoly](https://cran.r-project.org/package=fitPoly), and
[TET](http://www.g3journal.org/content/suppl/2017/01/19/g3.117.039008.DC1),
and [polyRAD](https://cran.r-project.org/package=polyRAD). Our best
“competitor” is probably
[fitPoly](https://cran.r-project.org/package=fitPoly), though
[polyRAD](https://cran.r-project.org/package=polyRAD) has some nice
ideas for utilizing population structure and linkage disequilibrium.

See [NEWS](https://github.com/dcgerard/updog/blob/master/NEWS.md) for
the latest updates on the package.

## Vignettes

I’ve included many vignettes in `updog`, which you can access online
[here](https://dcgerard.github.io/updog/).

## Bug Reports

If you find a bug or want an enhancement, please submit an issue
[here](http://github.com/dcgerard/updog/issues).

## Installation

You can install updog from
[CRAN](https://cran.r-project.org/package=updog) in the usual way:

``` r
install.packages("updog")
```

You can install the current (unstable) version of updog from GitHub
with:

``` r
# install.packages("devtools")
devtools::install_github("dcgerard/updog")
```

## How to Cite

Please cite

> Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018).
> Genotyping Polyploids from Messy Sequencing Data. *Genetics*, 210(3),
> 789-807. doi:
> [10.1534/genetics.118.301468](https://doi.org/10.1534/genetics.118.301468).

Or, using BibTex:

``` tex
@article {gerard2018genotyping,
    author = {Gerard, David and Ferr{\~a}o, Lu{\'i}s Felipe Ventorim and Garcia, Antonio Augusto Franco and Stephens, Matthew},
    title = {Genotyping Polyploids from Messy Sequencing Data},
    volume = {210},
    number = {3},
    pages = {789--807},
    year = {2018},
    doi = {10.1534/genetics.118.301468},
    publisher = {Genetics},
    issn = {0016-6731},
    URL = {https://doi.org/10.1534/genetics.118.301468},
    journal = {Genetics}
}
```

If you are using the proportional normal prior class (`model = "norm"`)
or the unimodal prior class (`model = "ash"`), then please also cite

> Gerard, D. & Ferrão L. F. V. (2019). “Priors for Genotyping
> Polyploids.” *Bioinformatics* (in press). doi:
> [10.1093/bioinformatics/btz852](https://doi.org/10.1093/bioinformatics/btz852)

Or, using BibTex:

``` tex
@article{gerard2019priors,
    author = {Gerard, David and Ferr{\~a}o, Lu{\'i}s Felipe Ventorim},
    title = {Priors for Genotyping Polyploids},
    journal = {Bioinformatics},
    year = {2019},
    month = {11},
    issn = {1367-4803},
    doi = {10.1093/bioinformatics/btz852},
    note = {btz852},
}
```

## Code of Conduct

Please note that this project is released with a [Contributor Code of
Conduct](https://github.com/dcgerard/updog/blob/master/CONDUCT.md). By
participating in this project you agree to abide by its terms.
