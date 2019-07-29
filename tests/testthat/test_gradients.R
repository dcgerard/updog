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

test_that("grad_for_mu_sigma2 doesn't return nan's", {
  set.seed(1)
  nind <- 11
  mu <- rnorm(nind)
  sigma2 <- abs(rnorm(nind))
  alpha <- 0.1
  rho   <- runif(nind)
  ploidy <- 6
  phifk_mat <- compute_all_phifk(alpha, rho, ploidy)[, 1, ]
  cor_inv <- solve(cov2cor(crossprod(matrix(rnorm(nind ^ 2), nrow = nind))))
  log_bb_dense <- matrix(abs(rnorm(nind * (ploidy + 1))), nrow = nind) * -100

  phifk_mat[, ploidy + 1] <- Inf
  phifk_mat[, 2] <- -Inf

  expect_false(any(is.nan(grad_for_mu_sigma2_wrapper(muSigma2 = c(mu, sigma2), phifk_mat = phifk_mat, cor_inv = cor_inv, log_bb_dense = log_bb_dense))))
}
)

test_that("dpen_dh works", {
  h        <- 1.3
  mu_h     <- 0.5
  sigma2_h <- 0.4

  dout <- dpen_dh(h = h, mu_h = mu_h, sigma2_h = sigma2_h)

  myenv <- new.env()
  assign(x = "h", value = h, envir = myenv)
  assign(x = "mu_h", value = mu_h, envir = myenv)
  assign(x = "sigma2_h", value = sigma2_h, envir = myenv)
  nout <- stats::numericDeriv(quote(pen_bias(h = h, mu_h = mu_h, sigma2_h = sigma2_h)), "h", myenv)
  expect_equal(attr(nout, "gradient")[1], dout)
}
)

test_that("dpen_deps works", {
  eps <- 0.01
  mu_eps <- -4.7
  sigma2_eps <- 1
  dout <- dpen_deps(eps = eps, mu_eps = mu_eps, sigma2_eps = sigma2_eps)

  myenv <- new.env()
  assign(x = "eps", value = eps, envir = myenv)
  assign(x = "mu_eps", value = mu_eps, envir = myenv)
  assign(x = "sigma2_eps", value = sigma2_eps, envir = myenv)
  nout <- stats::numericDeriv(quote(pen_seq_error(eps = eps, mu_eps = mu_eps, sigma2_eps = sigma2_eps)), "eps", myenv)

  expect_equal(attr(nout, "gradient")[1], dout)
}
)

test_that("dlbeta_dxi works", {
  x   <- 3
  n   <- 13
  xi  <- 0.3
  tau <- 0.1
  dout <- dlbeta_dxi(x = x, n = n, xi = xi, tau = tau)

  myenv <- new.env()
  assign(x = "x", value = x, envir = myenv)
  assign(x = "n", value = n, envir = myenv)
  assign(x = "xi", value = xi, envir = myenv)
  assign(x = "tau", value = tau, envir = myenv)
  nout <- stats::numericDeriv(quote(dbetabinom_double(x = x, size = n, mu = xi, rho = tau, log = TRUE)), "xi", myenv)

  expect_equal(attr(nout, "gradient")[1], dout, tol = 10^-6)
}
)

test_that("dlbeta_dh, dlbeta_deps, and dlbeta_dtau work", {
  x <- 3
  n <- 13
  p <- 1/6
  eps <- 0.02
  h <- 0.5
  tau <- 0.01

  dout1 <- dlbeta_dh(x = x, n = n, p = p, eps = eps, h = h, tau = tau)
  dout2 <- dlbeta_deps(x = x, n = n, p = p, eps = eps, h = h, tau = tau)
  dout3 <- dlbeta_dtau(x = x, n = n, p = p, eps = eps, h = h, tau = tau)

  fn <- function(x, n, p, eps, h, tau) {
    xi <- xi_double(p = p, eps = eps, h = h)
    dbetabinom_double(x = x, size = n, mu = xi, rho = tau, log = TRUE)
  }

  myenv <- new.env()
  assign(x = "x", value = x, envir = myenv)
  assign(x = "n", value = n, envir = myenv)
  assign(x = "p", value = p, envir = myenv)
  assign(x = "eps", value = eps, envir = myenv)
  assign(x = "h", value = h, envir = myenv)
  assign(x = "tau", value = tau, envir = myenv)
  nout1 <- stats::numericDeriv(quote(fn(x = x, n = n, p = p, eps = eps, h = h, tau = tau)), "h", myenv)
  nout2 <- stats::numericDeriv(quote(fn(x = x, n = n, p = p, eps = eps, h = h, tau = tau)), "eps", myenv)
  nout3 <- stats::numericDeriv(quote(fn(x = x, n = n, p = p, eps = eps, h = h, tau = tau)), "tau", myenv)

  expect_equal(attr(nout1, "gradient")[1], dout1, tol = 10 ^ -4)
  expect_equal(attr(nout2, "gradient")[1], dout2, tol = 10 ^ -4)
  expect_equal(attr(nout3, "gradient")[1], dout3, tol = 10 ^ -4)
}
)


test_that("grad_for_eps works", {
  set.seed(1)
  parvec  <- c(0.02, 1.5, 0.01)
  sizevec <- rpois(n = 11, lambda = 20)
  refvec  <- rbinom(n = 11, size = sizevec, prob = 0.4)
  ploidy <- 4
  mean_bias <- 0
  var_bias <- 1
  mean_seq <- -4.7
  var_seq <- 1
  mean_od <- -5.5
  var_od <- 0.6
  wmat <- matrix(abs(rnorm(11 * (ploidy + 1))), nrow = 11)
  wmat <- wmat / rowSums(wmat)

  refvec[1]  <- NA
  sizevec[2] <- NA

  dout <- grad_for_eps(parvec = parvec, refvec = refvec, sizevec = sizevec,
                       ploidy = ploidy, mean_bias = mean_bias,
                       var_bias = var_bias, mean_seq = mean_seq,
                       var_seq = var_seq, mean_od = mean_od,
                       var_od = var_od, wmat = wmat)

  myenv <- new.env()
  assign(x = "parvec", value = parvec, envir = myenv)
  assign(x = "refvec", value = refvec, envir = myenv)
  assign(x = "sizevec", value = sizevec, envir = myenv)
  assign(x = "ploidy", value = ploidy, envir = myenv)
  assign(x = "mean_bias", value = mean_bias, envir = myenv)
  assign(x = "var_bias", value = var_bias, envir = myenv)
  assign(x = "mean_seq", value = mean_seq, envir = myenv)
  assign(x = "var_seq", value = var_seq, envir = myenv)
  assign(x = "wmat", value = wmat, envir = myenv)
  nout <- stats::numericDeriv(quote(obj_for_eps(parvec = parvec, refvec = refvec, sizevec = sizevec,
                                                ploidy = ploidy, mean_bias = mean_bias,
                                                var_bias = var_bias, mean_seq = mean_seq,
                                                var_seq = var_seq, mean_od = mean_od,
                                                var_od = var_od, wmat = wmat)),
                              "parvec", myenv)
  expect_equal(c(attr(nout, "gradient")), dout, tol = 10 ^ -4)

  ## About twice as fast.
  # microbenchmark::microbenchmark(
  #   dout <- grad_for_eps(parvec = parvec, refvec = refvec, sizevec = sizevec,
  #                        ploidy = ploidy, mean_bias = mean_bias,
  #                        var_bias = var_bias, mean_seq = mean_seq,
  #                        var_seq = var_seq, wmat = wmat),
  #   nout <- stats::numericDeriv(quote(obj_for_eps(parvec = parvec, refvec = refvec, sizevec = sizevec,
  #                                                 ploidy = ploidy, mean_bias = mean_bias,
  #                                                 var_bias = var_bias, mean_seq = mean_seq,
  #                                                 var_seq = var_seq, wmat = wmat)),
  #                               "parvec", myenv)
  # )
}
)


test_that("grad_for_weighted_lbb is gradient for obj_for_weighted_lbb", {
  set.seed(1)
  ploidy <- 6
  weight_vec <- runif(ploidy + 1)
  parvec <- c(0.3, 0.1)

  myg <- grad_for_weighted_lbb(parvec = parvec, ploidy = ploidy, weight_vec = weight_vec)

  myenv <- new.env()
  assign(x = "parvec", value = parvec, envir = myenv)
  assign(x = "ploidy", value = ploidy, envir = myenv)
  assign(x = "weight_vec", value = weight_vec, envir = myenv)
  nout <- stats::numericDeriv(quote(obj_for_weighted_lbb(parvec = parvec,
                                                         ploidy = ploidy,
                                                         weight_vec = weight_vec)),
                              "parvec", myenv)
  expect_equal(c(attr(nout, "gradient")), myg, tol = 10^-5)
})


test_that("grad_for_weighted_lnorm is gradient for obj_for_weighted_lnorm", {
  set.seed(1)
  ploidy <- 6
  weight_vec <- runif(ploidy + 1)
  parvec <- c(11, 2)

  myg <- grad_for_weighted_lnorm(parvec = parvec, ploidy = ploidy, weight_vec = weight_vec)

  myenv <- new.env()
  assign(x = "parvec", value = parvec, envir = myenv)
  assign(x = "ploidy", value = ploidy, envir = myenv)
  assign(x = "weight_vec", value = weight_vec, envir = myenv)
  nout <- stats::numericDeriv(quote(obj_for_weighted_lnorm(parvec = parvec,
                                                           ploidy = ploidy,
                                                           weight_vec = weight_vec)),
                              "parvec", myenv)
  expect_equal(c(attr(nout, "gradient")), myg, tol = 10^-5)


})




