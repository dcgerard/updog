# EM algorithm to fit weighted ash objective.

Solves the following optimization problem \$\$\max\_{\pi} \sum_k w_k
\log(\sum_j \pi_j \ell_jk).\$\$ It does this using a weighted EM
algorithm.

## Usage

``` r
wem(weight_vec, lmat, pi_init, lambda, itermax, obj_tol)
```

## Arguments

- weight_vec:

  A vector of weights. Each element of `weight_vec` corresponds to a
  column of `lmat`.

- lmat:

  A matrix of inner weights. The columns are the "individuals" and the
  rows are the "classes."

- pi_init:

  The initial values of `pivec`. Each element of `pi_init` corresponds
  to a row of `lmat`.

- lambda:

  The penalty on the pi's. Should be greater than 0 and really really
  small.

- itermax:

  The maximum number of EM iterations to take.

- obj_tol:

  The objective stopping criterion.

## Value

A vector of numerics.

## Author

David Gerard

## Examples

``` r
set.seed(2)
n <- 3
p <- 5
lmat <- matrix(stats::runif(n * p), nrow = n)
weight_vec <- seq_len(p)
pi_init <- stats::runif(n)
pi_init <- pi_init / sum(pi_init)
wem(weight_vec = weight_vec,
    lmat       = lmat,
    pi_init    = pi_init,
    lambda     = 0,
    itermax    = 100,
    obj_tol    = 10^-6)
#>              [,1]
#> [1,] 3.830930e-01
#> [2,] 6.169070e-01
#> [3,] 3.041614e-09

```
