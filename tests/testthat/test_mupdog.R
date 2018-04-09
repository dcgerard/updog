context("mupdog")

test_that("mupdog works", {
  set.seed(1)
  nind <- 11
  nsnps <- 37

  sizemat <- matrix(rpois(n = nind * nsnps, lambda = 100),
                    nrow = nind, ncol = nsnps)
  refmat <- matrix(rbinom(n = nind * nsnps,
                          size = c(sizemat),
                          prob = rep(seq(0.1, 0.8, length = nsnps),
                                     each = nind)),
                   nrow = nind, ncol = nsnps)

  seq                <- rep(0.005, length = nsnps)
  bias               <- rep(1, length = nsnps)
  od                 <- rep(0.001, length = nsnps)
  allele_freq        <- colMeans(refmat / sizemat, na.rm = TRUE)
  inbreeding         <- rep(0.001, length = nind)
  cor_mat            <- diag(nind)
  postmean           <- matrix(0, nrow = nind, ncol = nsnps)
  postvar            <- matrix(1, nrow = nind, ncol = nsnps)
  var_bias           <- 1
  mean_bias          <- 0
  var_seq            <- 1
  mean_seq           <- -4.7
  update_cor         <- TRUE
  update_inbreeding  <- TRUE
  update_allele_freq <- TRUE
  control            <- list()
  ploidy             <- 4
  verbose            <- FALSE
  num_core           <- 2
  itermax            <- 100
  obj_tol            <- 10 ^ -4

  refmat[1,1] <- NA

  trash <- capture.output(
    mout <- mupdog(refmat = refmat, sizemat = sizemat, ploidy = 8,
                   mean_bias = 0, var_bias = 1, mean_seq = -4.7, var_seq = 1,
                   seq = seq, bias = bias, od = od,
                   allele_freq = allele_freq,
                   inbreeding = inbreeding,
                   cor_mat = cor_mat,
                   postmean = postmean,
                   postvar = postvar,
                   control = list(itermax = 2), verbose = FALSE)
  )

  warray <- compute_all_post_prob(ploidy = 8, mu = postmean,
                                  sigma2 = postvar,
                                  alpha = allele_freq,
                                  rho = inbreeding)
  expect_equal(sum(warray[10,2,]), 1)
  expect_true(all(warray > 0))
  expect_true(!is.nan(post_prob(8, 8, 0, 1, 0.1, 0.001)))
  expect_true(pbetabinom_double(8, 8, 0.1, 0.001, FALSE) <= 1)

}
)
