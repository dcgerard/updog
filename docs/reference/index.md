# Package index

## All functions

- [`dbetabinom()`](https://dcgerard.github.io/updog/reference/betabinom.md)
  [`pbetabinom()`](https://dcgerard.github.io/updog/reference/betabinom.md)
  [`qbetabinom()`](https://dcgerard.github.io/updog/reference/betabinom.md)
  [`rbetabinom()`](https://dcgerard.github.io/updog/reference/betabinom.md)
  : The Beta-Binomial Distribution

- [`filter_snp()`](https://dcgerard.github.io/updog/reference/filter_snp.md)
  :

  Filter SNPs based on the output of
  [`multidog()`](https://dcgerard.github.io/updog/reference/multidog.md).

- [`flexdog()`](https://dcgerard.github.io/updog/reference/flexdog.md) :
  Flexible genotyping for polyploids from next-generation sequencing
  data.

- [`flexdog_full()`](https://dcgerard.github.io/updog/reference/flexdog_full.md)
  : Flexible genotyping for polyploids from next-generation sequencing
  data.

- [`format_multidog()`](https://dcgerard.github.io/updog/reference/format_multidog.md)
  :

  Return arrayicized elements from the output of `multidog`.

- [`get_q_array()`](https://dcgerard.github.io/updog/reference/get_q_array.md)
  : Return the probabilities of an offspring's genotype given its
  parental genotypes for all possible combinations of parental and
  offspring genotypes. This is for species with polysomal inheritance
  and bivalent, non-preferential pairing.

- [`is.flexdog()`](https://dcgerard.github.io/updog/reference/is.flexdog.md)
  :

  Tests if an argument is a `flexdog` object.

- [`is.multidog()`](https://dcgerard.github.io/updog/reference/is.multidog.md)
  :

  Tests if an argument is a `multidog` object.

- [`log_sum_exp()`](https://dcgerard.github.io/updog/reference/log_sum_exp.md)
  : Log-sum-exponential trick.

- [`log_sum_exp_2()`](https://dcgerard.github.io/updog/reference/log_sum_exp_2.md)
  : Log-sum-exponential trick using just two doubles.

- [`multidog()`](https://dcgerard.github.io/updog/reference/multidog.md)
  :

  Fit `flexdog` to multiple SNPs.

- [`oracle_cor()`](https://dcgerard.github.io/updog/reference/oracle_cor.md)
  : Calculates the correlation between the true genotype and an oracle
  estimator.

- [`oracle_cor_from_joint()`](https://dcgerard.github.io/updog/reference/oracle_cor_from_joint.md)
  : Calculate the correlation of the oracle estimator with the true
  genotype from the joint distribution matrix.

- [`oracle_joint()`](https://dcgerard.github.io/updog/reference/oracle_joint.md)
  : The joint probability of the genotype and the genotype estimate of
  an oracle estimator.

- [`oracle_mis()`](https://dcgerard.github.io/updog/reference/oracle_mis.md)
  : Calculate oracle misclassification error rate.

- [`oracle_mis_from_joint()`](https://dcgerard.github.io/updog/reference/oracle_mis_from_joint.md)
  : Get the oracle misclassification error rate directly from the joint
  distribution of the genotype and the oracle estimator.

- [`oracle_mis_vec()`](https://dcgerard.github.io/updog/reference/oracle_mis_vec.md)
  : Returns the oracle misclassification rates for each genotype.

- [`oracle_mis_vec_from_joint()`](https://dcgerard.github.io/updog/reference/oracle_mis_vec_from_joint.md)
  : Get the oracle misclassification error rates (conditional on true
  genotype) directly from the joint distribution of the genotype and the
  oracle estimator.

- [`oracle_plot()`](https://dcgerard.github.io/updog/reference/oracle_plot.md)
  :

  Construct an oracle plot from the output of `oracle_joint`.

- [`plot(`*`<flexdog>`*`)`](https://dcgerard.github.io/updog/reference/plot.flexdog.md)
  :

  Draw a genotype plot from the output of `flexdog`.

- [`plot(`*`<multidog>`*`)`](https://dcgerard.github.io/updog/reference/plot.multidog.md)
  :

  Plot the output of `multidog`.

- [`plot_geno()`](https://dcgerard.github.io/updog/reference/plot_geno.md)
  : Make a genotype plot.

- [`rflexdog()`](https://dcgerard.github.io/updog/reference/rflexdog.md)
  :

  Simulate GBS data from the `flexdog` likelihood.

- [`rgeno()`](https://dcgerard.github.io/updog/reference/rgeno.md) :

  Simulate individual genotypes from one of the supported `flexdog`
  models.

- [`snpdat`](https://dcgerard.github.io/updog/reference/snpdat.md) : GBS
  data from Shirasawa et al (2017)

- [`uitdewilligen`](https://dcgerard.github.io/updog/reference/uitdewilligen.md)
  : Subset of individuals and SNPs from Uitdewilligen et al (2013).

- [`wem()`](https://dcgerard.github.io/updog/reference/wem.md) : EM
  algorithm to fit weighted ash objective.
