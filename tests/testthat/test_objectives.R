context("objectives")

test_that("post_prob works", {

  cppway <- post_prob(2, 6, -1, 1, 0.5, 0.01)
  rway <- pnorm(qnorm(pbetabinom(2, 6, 0.5, 0.01, FALSE)), mean = -1, sd = 1) -
    pnorm(qnorm(pbetabinom(1, 6, 0.5, 0.01, FALSE)), mean = -1, sd = 1)

  expect_equal(cppway, rway)
}
)

test_that("pen_bias works", {
  cppway <- pen_bias(h = 0.8, mu_h = 1, sigma2_h = 1)
  rway   <- -log(0.8) - (log(0.8) - 1) ^ 2 / (2 * 1)
  expect_equal(cppway, rway)
}
)

test_that("pen_seq_error works", {
  cppway <- pen_seq_error(0.2, -4, 1)
  rway   <- -log(0.2 * 0.8) - (log(0.2 / 0.8) + 4) ^ 2 / (2 * 1)
  expect_equal(cppway, rway)
}
)

test_that("compute_all_log_bb returns NA", {

  refmat <- matrix(c(1,2,3,NA), nrow = 2)
  sizemat <- matrix(c(NA,3,4,5), nrow = 2)
  ploidy <- 4
  seq <- c(0.2, 0.1)
  bias <- c(1, 1.3)
  od <- c(0.01, 0.02)

  bbdenseout <- compute_all_log_bb(refmat = refmat, sizemat = sizemat, ploidy = ploidy, seq = seq, bias = bias, od = od)

  expect_true(all(is.na(bbdenseout[1, 1, ])))
  expect_true(all(is.na(bbdenseout[2, 2, ])))
  expect_false(all(is.na(bbdenseout[1, 2, ])))
  expect_false(all(is.na(bbdenseout[2, 1, ])))
}
)

test_that("compute_all_phifk works", {
  alpha <- c(0.1, 0.2)
  rho   <- c(0.01, 0.02, 0.03)
  ploidy <- 4
  phifk <- compute_all_phifk(alpha, rho, ploidy)[1, 1, , drop = TRUE]
  rway <- qnorm(pbetabinom(-1:ploidy, size = rep(ploidy, ploidy + 2), rho = rho[1], mu = alpha[1], log_p = FALSE))
  expect_equal(phifk[2:(ploidy + 1)], rway[2:(ploidy + 1)])
}
)

test_that("obj_for_mu_sigma2 and elbo match", {
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

  seq <- rep(0.005, length = nsnps)
  bias <- rep(1, length = nsnps)
  od <- rep(0.001, length = nsnps)
  allele_freq <- colMeans(refmat / sizemat, na.rm = TRUE)
  inbreeding <- rep(0.001, length = nind)
  cor_inv <- solve(cov2cor(crossprod(matrix(rnorm(nind ^ 2), nrow = nind))))
  postmean <- matrix(0, nrow = nind, ncol = nsnps)
  postvar <- matrix(1, nrow = nind, ncol = nsnps)
  var_bias = 1
  mean_bias = 0
  var_seq = 1
  mean_seq = -4.7
  ploidy <- 6

  warray <- compute_all_post_prob(ploidy = ploidy, mu = postmean, sigma2 = postvar,
                                  alpha = allele_freq, rho = inbreeding)
  lbeta_array <- compute_all_log_bb(refmat = refmat, sizemat = sizemat,
                                    ploidy = ploidy, seq = seq, bias = bias, od = od)
  phifk_array <- compute_all_phifk(alpha = allele_freq, rho = inbreeding, ploidy = ploidy)
  obj_elbo <- elbo(warray = warray, lbeta_array = lbeta_array, cor_inv = cor_inv,
                   postmean = postmean, postvar = postvar, bias = bias, seq = seq,
                   mean_bias = mean_bias, var_bias = var_bias, mean_seq = mean_seq,
                   var_seq = var_seq, ploidy = ploidy)

  obj_manual <- 0
  for (index in 1:nsnps) {
    obj_manual <- obj_manual + obj_for_mu_sigma2(mu = postmean[, index], sigma2 = postvar[, index], phifk_mat = phifk_array[, index, ],
                                                 cor_inv = cor_inv, log_bb_dense = lbeta_array[, index, ])
  }

  obj_manual <- obj_manual - sum((log(bias) - mean_bias) ^ 2) / (2 * var_bias) -
    sum(log(bias)) -
    sum(log(seq * (1 - seq))) -
    sum((log(seq / (1 - seq)) - mean_seq) ^ 2 / (2 * var_seq)) +
    determinant(cor_inv, logarithm = TRUE)$modulus[1] * nsnps / 2
  expect_equal(obj_manual, obj_elbo)


  ## test two ways to calculate obj_for_mu_sigma2
  obj1 <- obj_for_mu_sigma2(mu = postmean[, index], sigma2 = postvar[, index], phifk_mat = phifk_array[, index, ],
                            cor_inv = cor_inv, log_bb_dense = lbeta_array[, index, ])
  obj2 <- obj_for_mu_sigma2_wrapper(muSigma2 = c(postmean[, index], postvar[, index]), phifk_mat = phifk_array[, index, ],
                            cor_inv = cor_inv, log_bb_dense = lbeta_array[, index, ])
  expect_equal(obj1, obj2)
}
)

test_that("obj_for_weighted_lnorm is correct", {
  weight_vec <- c(1, 2, 3)
  ploidy <- 2
  parvec <- c(0.3, 0.4)

  cobj <- obj_for_weighted_lnorm(parvec = parvec, ploidy = ploidy, weight_vec = weight_vec)

  pivec <- dnorm(x = 0:ploidy, mean = parvec[1], sd = parvec[2])
  pivec <- pivec / sum(pivec)
  robj <- sum(log(pivec) * weight_vec)

  expect_equal(cobj, robj)

  ## The R code is actually a little faster
  # microbenchmark::microbenchmark(
  #   {pivec <- dnorm(x = 0:ploidy, mean = parvec[1], sd = parvec[2])
  #   pivec <- pivec / sum(pivec)
  #   robj <- sum(log(pivec) * weight_vec)},
  #   cobj <- obj_for_weighted_lnorm(parvec = parvec, ploidy = ploidy, weight_vec = weight_vec)
  # )
})

test_that("no penalty when variance is Inf", {
  expect_equal(pen_bias(h = 0.1, mu_h = 1, sigma2_h = Inf), 0)
  expect_equal(pen_seq_error(eps = 0.1, mu_eps = -1, sigma2_eps = Inf), 0)
  expect_equal(dpen_deps(eps = 0.1, mu_eps = -1, sigma2_eps = Inf), 0)  
  expect_equal(dpen_dh(h = 0.1, mu_h = 1, sigma2_h = Inf), 0)  
  
  expect_true(pen_bias(h = 0.1, mu_h = 1, sigma2_h = 10000) != 0)
  expect_true(pen_seq_error(eps = 0.1, mu_eps = -1, sigma2_eps = 10000) != 0)
  expect_true(dpen_deps(eps = 0.1, mu_eps = -1, sigma2_eps = 10000) != 0)  
  expect_true(dpen_dh(h = 0.1, mu_h = 1, sigma2_h = 10000) != 0)  
  
})
