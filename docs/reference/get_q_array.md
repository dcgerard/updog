# Return the probabilities of an offspring's genotype given its parental genotypes for all possible combinations of parental and offspring genotypes. This is for species with polysomal inheritance and bivalent, non-preferential pairing.

Return the probabilities of an offspring's genotype given its parental
genotypes for all possible combinations of parental and offspring
genotypes. This is for species with polysomal inheritance and bivalent,
non-preferential pairing.

## Usage

``` r
get_q_array(ploidy)
```

## Arguments

- ploidy:

  A positive integer. The ploidy of the species.

## Value

An three-way array of proportions. The (i, j, k)th element is the
probability of an offspring having k - 1 reference alleles given that
parent 1 has i - 1 reference alleles and parent 2 has j - 1 reference
alleles. Each dimension of the array is `ploidy + 1`. In the dimension
names, "A" stands for the reference allele and "a" stands for the
alternative allele.

## Author

David Gerard

## Examples

``` r
qarray <- get_q_array(6)
apply(qarray, c(1, 2), sum) ## should all be 1's.
#>         parent2
#> parent1  aaaaaa Aaaaaa AAaaaa AAAaaa AAAAaa AAAAAa AAAAAA
#>   aaaaaa      1      1      1      1      1      1      1
#>   Aaaaaa      1      1      1      1      1      1      1
#>   AAaaaa      1      1      1      1      1      1      1
#>   AAAaaa      1      1      1      1      1      1      1
#>   AAAAaa      1      1      1      1      1      1      1
#>   AAAAAa      1      1      1      1      1      1      1
#>   AAAAAA      1      1      1      1      1      1      1
```
