# Return arrayicized elements from the output of [`multidog`](https://dcgerard.github.io/updog/reference/multidog.md).

This function will allow you to have genotype estimates, maximum
posterior probability, and other values in the form of a matrix/array.
If multiple variable names are provided, the data are formatted as a
3-dimensional array with the dimensions corresponding to (individuals,
SNPs, variables).

## Usage

``` r
format_multidog(x, varname = "geno")
```

## Arguments

- x:

  The output of `multidog`.

- varname:

  A character vector of the variable names whose values populate the
  cells. These should be column names from `x$inddf`.

## Details

Note that the order of the individuals will be reshuffled. The order of
the SNPs should be the same as in `x$snpdf`.

## Author

David Gerard
