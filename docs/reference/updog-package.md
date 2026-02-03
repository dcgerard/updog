# updog: Flexible Genotyping for Polyploids

Implements empirical Bayes approaches to genotype polyploids from next
generation sequencing data while accounting for allele bias,
overdispersion, and sequencing error. The main functions are flexdog()
and multidog(), which allow the specification of many different genotype
distributions. Also provided are functions to simulate genotypes,
rgeno(), and read-counts, rflexdog(), as well as functions to calculate
oracle genotyping error rates, oracle_mis(), and correlation with the
true genotypes, oracle_cor(). These latter two functions are useful for
read depth calculations. Run browseVignettes(package = "updog") in R for
example usage. See Gerard et al. (2018)
[doi:10.1534/genetics.118.301468](https://doi.org/10.1534/genetics.118.301468)
and Gerard and Ferrao (2020)
[doi:10.1093/bioinformatics/btz852](https://doi.org/10.1093/bioinformatics/btz852)
for details on the implemented methods.

## See also

Useful links:

- <https://dcgerard.github.io/updog/>

- Report bugs at <https://github.com/dcgerard/updog/issues>

## Author

**Maintainer**: David Gerard <gerard.1787@gmail.com>
([ORCID](https://orcid.org/0000-0001-9450-5023))
