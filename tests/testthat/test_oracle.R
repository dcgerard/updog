context("oracle")

test_that("oracle error rates are correct", {
 n      <- 100
 ploidy <- 6
 seq    <- 0.01
 bias   <- 0.75
 od     <- 0
 dist   <- dbinom(x = 0:ploidy, size = ploidy, prob = 0.7)
 oracle_mis(n = n, ploidy = ploidy, seq = seq, bias = bias, od = od, dist = dist)


 rfun <- function(n, ploidy, seq, bias, od, dist) {
   xi <- xi_fun(p = (0:ploidy) / ploidy, eps = seq, h = bias)

   dbbmat <- matrix(NA, nrow = n + 1, ncol = ploidy + 1)
   for (index in 1:(ploidy + 1)) {
     dbbmat[, index] <- dist[index] * dbetabinom(x = 0:n, size = n, mu = xi[index], rho = od, log = FALSE)
   }

   1 - sum(apply(dbbmat, 1, max))

 }

 rway <- rfun(n = n, ploidy = ploidy, seq = seq, bias = bias, od = od, dist = dist)
 myway <- oracle_mis(n = n, ploidy = ploidy, seq = seq, bias = bias, od = od, dist = dist)

 expect_equal(rway, myway)

 # ## 7 times faster!
 # microbenchmark::microbenchmark(
 #   rfun(n = n, ploidy = ploidy, seq = seq, bias = bias, od = od, dist = dist),
 #   oracle_mis(n = n, ploidy = ploidy, seq = seq, bias = bias, od = od, dist = dist)
 # )

}
)
