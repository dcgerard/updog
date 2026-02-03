# Filter SNPs based on the output of [`multidog()`](https://dcgerard.github.io/updog/reference/multidog.md).

Filter based on provided logical predicates in terms of the variable
names in `x$snpdf`. This function filters both `x$snpdf` and `x$inddf`.

## Usage

``` r
filter_snp(x, expr)
```

## Arguments

- x:

  The output of `multidog`.

- expr:

  Logical predicate expression defined in terms of the variables in
  `x$snpdf`. Only SNPs where the condition evaluates to `TRUE` are kept.

## See also

- [`multidog()`](https://dcgerard.github.io/updog/reference/multidog.md)::

  For the variables in `x$snpdf` which you can filter by.

## Author

David Gerard

## Examples

``` r
# \donttest{
data("uitdewilligen")
mout <- multidog(refmat = t(uitdewilligen$refmat),
                 sizemat = t(uitdewilligen$sizemat),
                 ploidy = uitdewilligen$ploidy,
                 nc = 2)
#>     |                                   *.#,%    
#>    |||                                 *******/  
#>  |||||||    (**..#**.                  */   **/  
#> |||||||||    */****************************/*%   
#>    |||    &****..,*.************************/    
#>    |||     (....,,,*,...****%********/(******    
#>    |||                ,,****%////,,,,./.****/    
#>    |||                  /**//         .*///....  
#>    |||                  .*/*/%#         .,/   ., 
#>    |||               , **/   #%         .*    .. 
#>    |||                               ,,,*        
#> 
#> Working on it...done!

## The following filters are for educational purposes only and should
## not be taken as a default filter:
mout2 <- filter_snp(mout, bias < 0.8 & od < 0.003)
# }
```
