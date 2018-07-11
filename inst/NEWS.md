# Errata

* The original documentation of `flexdog` says

    > The value of `prop_mis` is a very intuitive measure for the quality of the SNP. `prop_mis` is the posterior proportion of individuals mis-genotyped. So if you want only SNPS that accurately genotype, say, 95% of the individuals, you could discard all SNPs with a `prop_mis` under 0.95.

    This now says
    
    > The value of prop_mis is a very intuitive measure for the quality of the SNP. prop_mis is the posterior proportion of individuals mis-genotyped. So if you want only SNPS that accurately genotype, say, 95% of the individuals, you could discard all SNPs with a prop_mis **over 0.05**.

* The original documentation of `flexdog` also says
  
    > The value of `maxpostprob` is a very intuitive measure for the quality of the genotype estimate of an individual. This is the posterior probability of correctly genotyping the individual when using `geno` (the posterior mode) as the genotype estimate. So if you want to correctly genotype, say, 95% of individuals, you could discard all individuals with a `maxpostprob` of under 0.95.
    
    This should now says
    
    > > The value of `maxpostprob` is a very intuitive measure for the quality of the genotype estimate of an individual. This is the posterior probability of correctly genotyping the individual when using `geno` (the posterior mode) as the genotype estimate. So if you want to correctly genotype, say, 95% of individuals, you could discard all individuals with a `maxpostprob` of **over 0.05**.


# updog 0.99.0

* This is a complete re-working of the code in `updog`. The old version may be found in the [`updogAlpha` package](https://github.com/dcgerard/updogAlpha).
* The main function is now `flexdog`.
* An experimental approach  `mupdog` is now live. We provide no guarantees about `mupdog`'s performance.
* Oracle misclassification error rates may be calculated in `oracle_mis`.
* Genotypes can be simulated using `rgeno`.
* Next-generation sequencing data can be simulated using `rflexdog`.
