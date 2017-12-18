context("gradients")

test_that("grad_for_mu_sigma2 is the gradient for obj_for_mu_sigma2", {
  set.seed(1)
  nind <- 11
  mu <- rnorm(nind)
  sigma2 <- abs(rnorm(nind))
  alpha <- 0.1
  rho   <- runif(nind)
  ploidy <- 6
  phifk_mat <- compute_all_phifk(alpha, rho, ploidy)[, 1, ]
  cor_inv <- solve(cov2cor(crossprod(matrix(rnorm(nind ^ 2), nrow = nind))))
  # cor_inv <- diag(nind)
  log_bb_dense <- matrix(abs(rnorm(nind * (ploidy + 1))), nrow = nind) * -100

  obj_for_mu_sigma2(mu = mu, sigma2 = sigma2, phifk_mat = phifk_mat, cor_inv = cor_inv, log_bb_dense = log_bb_dense)
  gradvec <- grad_for_mu_sigma2(mu = mu, sigma2 = sigma2, phifk_mat = phifk_mat, cor_inv = cor_inv, log_bb_dense = log_bb_dense)

  myenv <- new.env()
  assign(x = "mu", value = mu, envir = myenv)
  assign(x = "sigma2", value = sigma2, envir = myenv)
  assign(x = "phifk_mat", value = phifk_mat, envir = myenv)
  assign(x = "cor_inv", value = cor_inv, envir = myenv)
  assign(x = "log_bb_dense", value = log_bb_dense, envir = myenv)
  nout <- stats::numericDeriv(quote(obj_for_mu_sigma2(mu, sigma2, phifk_mat, cor_inv, log_bb_dense)), "mu", myenv)
  expect_equal(c(attr(nout, "gradient")), gradvec[1:nind], tol = 10^-5)

  nout <- stats::numericDeriv(quote(obj_for_mu_sigma2(mu, sigma2, phifk_mat, cor_inv, log_bb_dense)), "sigma2", myenv)
  expect_equal(c(attr(nout, "gradient")), gradvec[(nind + 1): (2 * nind)], tol = 10^-5)


  ## now test if grad_for_mu_sigma2 is the same as grad_for_mu_sigma2_wrapper
  gradvec2 <- grad_for_mu_sigma2_wrapper(muSigma2 = c(mu, sigma2), phifk_mat = phifk_mat, cor_inv = cor_inv, log_bb_dense = log_bb_dense)
  expect_equal(gradvec, gradvec2)
}
)
