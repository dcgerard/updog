context("uni_em vs uni_gem")


test_that("uni_gem works", {
  n <- 7
  k <- 7
  lmat <- matrix(runif(n * k), nrow = k)
  weight_vec <- rep(1, length = k)
  pi_init <- rep(1 / k, k)

  piuni <- uni_em(weight_vec = weight_vec,
                  lmat       = lmat,
                  pi_init    = pi_init,
                  lambda     = 10^-8,
                  itermax    = 1000,
                  obj_tol    = 10^-5)
  piw <- wem(weight_vec = weight_vec,
             lmat       = lmat,
             pi_init    = pi_init,
             lambda     = 10^-8,
             itermax    = 1000,
             obj_tol    = 10^-5)

  trash <- wem_obj(pivec = piw,
                   weight_vec = weight_vec,
                   lmat = lmat,
                   lambda = 10^-8)
  trash <- wem_obj(pivec = piuni,
                   weight_vec = weight_vec,
                   lmat = lmat,
                   lambda = 10^-8)
})


test_that("uni_em_const gives same results as uni_em when alpha = 0", {
  set.seed(3)
  k <- 7
  n <- 7
  lmat <- matrix(runif(n * k), nrow = k)
  weight_vec <- rep(1, length = k)
  pi_init <- rep(1 / k, k)

  u1 <- uni_em(weight_vec = weight_vec,
               lmat       = lmat,
               pi_init    = pi_init,
               lambda     = 10^-8,
               itermax    = 1000,
               obj_tol    = 10^-5)

  u2 <- uni_em_const(weight_vec = weight_vec,
                     lmat       = lmat,
                     pi_init    = pi_init,
                     alpha      = 0,
                     lambda     = 10^-8,
                     itermax    = 1000,
                     obj_tol    = 10^-5)
  
  u3 <- uni_em_const(weight_vec = weight_vec,
                     lmat       = lmat,
                     pi_init    = pi_init,
                     alpha      = 0,
                     lambda     = c(rep(100, 3), rep(1, 4)),
                     itermax    = 1000,
                     obj_tol    = 10^-5)
  
  expect_error(uni_em_const(weight_vec = weight_vec,
                            lmat       = lmat,
                            pi_init    = pi_init,
                            alpha      = 0,
                            lambda     = c(rep(10^-8, 3), rep(1, )),
                            itermax    = 1000,
                            obj_tol    = 10^-5))

  expect_equal(u1, u2)

  ## make sure no errors
  u2 <- uni_em_const(weight_vec = weight_vec,
                     lmat       = lmat,
                     pi_init    = pi_init,
                     alpha      = 0.5,
                     lambda     = 10^-8,
                     itermax    = 1000,
                     obj_tol    = 10^-5)


  expect_equal(sum(u2), 1)


  k <- 3
  n <- 7
  lmat <- matrix(runif(n * k), nrow = k)
  weight_vec <- runif(n)
  pi_init <- rep(1 / k, k)
  u2 <- uni_em_const(weight_vec = weight_vec,
                     lmat       = lmat,
                     pi_init    = pi_init,
                     alpha      = 0,
                     lambda     = 10^-8,
                     itermax    = 1000,
                     obj_tol    = 10^-7)

  piw <- wem(weight_vec = weight_vec,
             lmat       = lmat,
             pi_init    = pi_init,
             lambda     = 10^-8,
             itermax    = 1000,
             obj_tol    = 10^-7)

  expect_equal(piw, u2, tol = 10^-3)


})
