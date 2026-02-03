# The Beta-Binomial Distribution

Density, distribution function, quantile function and random generation
for the beta-binomial distribution when parameterized by the mean `mu`
and the overdispersion parameter `rho` rather than the typical shape
parameters.

## Usage

``` r
dbetabinom(x, size, mu, rho, log)

pbetabinom(q, size, mu, rho, log_p)

qbetabinom(p, size, mu, rho)

rbetabinom(n, size, mu, rho)
```

## Arguments

- x, q:

  A vector of quantiles.

- size:

  A vector of sizes.

- mu:

  Either a scalar of the mean for each observation, or a vector of means
  of each observation, and thus the same length as `x` and `size`. This
  must be between 0 and 1.

- rho:

  Either a scalar of the overdispersion parameter for each observation,
  or a vector of overdispersion parameters of each observation, and thus
  the same length as `x` and `size`. This must be between 0 and 1.

- log, log_p:

  A logical vector either of length 1 or the same length as `x` and
  `size`. This determines whether to return the log probabilities for
  all observations (in the case that its length is 1) or for each
  observation (in the case that its length is that of `x` and `size`).

- p:

  A vector of probabilities.

- n:

  The number of observations.

## Value

Either a random sample (`rbetabinom`), the density (`dbetabinom`), the
tail probability (`pbetabinom`), or the quantile (`qbetabinom`) of the
beta-binomial distribution.

## Details

Let \\\mu\\ and \\\rho\\ be the mean and overdispersion parameters. Let
\\\alpha\\ and \\\beta\\ be the usual shape parameters of a beta
distribution. Then we have the relation \$\$\mu = \alpha/(\alpha +
\beta),\$\$ and \$\$\rho = 1/(1 + \alpha + \beta).\$\$ This necessarily
means that \$\$\alpha = \mu (1 - \rho)/\rho,\$\$ and \$\$\beta = (1 -
\mu) (1 - \rho)/\rho.\$\$

## Functions

- `dbetabinom()`: Density function.

- `pbetabinom()`: Distribution function.

- `qbetabinom()`: Quantile function.

- `rbetabinom()`: Random generation.

## Author

David Gerard

## Examples

``` r
x <- rbetabinom(n = 10, size = 10, mu = 0.1, rho = 0.01)
dbetabinom(x = 1, size = 10, mu = 0.1, rho = 0.01, log = FALSE)
#> [1] 0.3689335
pbetabinom(q = 1, size = 10, mu = 0.1, rho = 0.01, log_p = FALSE)
#> [1] 0.7345131
qbetabinom(p = 0.6, size = 10, mu = 0.1, rho = 0.01)
#> [1] 1

```
