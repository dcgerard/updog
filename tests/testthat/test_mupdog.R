context("mupdog")

test_that("mupdog works", {
  nind <- 11
  nsnps <- 37

  sizemat <- matrix(rpois(n = nind * nsnps, lambda = 100),
                    nrow = nind, ncol = nsnps)
  refmat <- matrix(rbinom(n = nind * nsnps,
                          size = c(sizemat),
                          prob = rep(seq(0.1, 0.8, length = nsnps),
                                     each = nind)),
                   nrow = nind, ncol = nsnps)

  seq_init <- rep(0.005, length = nsnps)
  bias_init <- rep(0, length = nsnps)
  od_init <- rep(0.001, length = nsnps)
  allele_freq_init <- colMeans(refmat / sizemat, na.rm = TRUE)
  inbreeding_init <- rep(0.001, length = nind)
  cor_init <- diag(nind)
  postmean_init <- matrix(0, nrow = nind, ncol = nsnps)
  postvar_init <- matrix(1, nrow = nind, ncol = nsnps)

  mupdog(refmat = refmat, sizemat = sizemat, ploidy = 8,
         mean_bias = 0, var_bias = 1, mean_seq = -4.7, var_seq = 1,
         seq_init = seq_init, bias_init = bias_init, od_init = od_init,
         allele_freq_init = allele_freq_init,
         inbreeding_init = inbreeding_init,
         cor_init = cor_init,
         postmean_init = postmean_init,
         postvar_init = postvar_init)

}
)
