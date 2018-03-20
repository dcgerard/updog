
<!-- README.md is generated from README.Rmd. Please edit that file -->
mupdog
======

[![Travis-CI Build Status](https://travis-ci.org/dcgerard/mupdog.svg?branch=master)](https://travis-ci.org/dcgerard/mupdog) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/dcgerard/mupdog?branch=master&svg=true)](https://ci.appveyor.com/project/dcgerard/mupdog) [![Coverage Status](https://img.shields.io/codecov/c/github/dcgerard/mupdog/master.svg)](https://codecov.io/github/dcgerard/mupdog?branch=master) [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Mupdog (multi-SNP updog) is an extension of the updog procedure for genotyping autopolyploids using GBS or RAD-seq like data. It allows for correlation between the individuals' genotypes and jointly estimates the genotypes of the individuals at all provided SNPs. The implementation uses a variational approximation.

See also [updog](https://github.com/dcgerard/updog), [ebg](https://github.com/pblischak/polyploid-genotyping), and [TET](http://www.g3journal.org/content/suppl/2017/01/19/g3.117.039008.DC1).

This package also contains `flexdog`, a function that flexibly estimates the genotype distribution.

Installation
------------

You can install mupdog from Github with:

``` r
# install.packages("devtools")
devtools::install_github("dcgerard/mupdog")
```

`flexdog` depends on the [CVXR](https://cran.r-project.org/web/packages/CVXR/index.html) package. Before I could install CVXR in Ubuntu, I had to run in the terminal

``` bash
sudo apt-get install libmpfr-dev
```

and then run in R

``` r
install.packages("Rmpfr")
```

Code of Conduct
---------------

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.
