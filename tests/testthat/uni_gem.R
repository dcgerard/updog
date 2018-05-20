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
  
  wem_obj(pivec = piw, 
          weight_vec = weight_vec, 
          lmat = lmat, 
          lambda = 10^-8)
  wem_obj(pivec = piuni, 
          weight_vec = weight_vec, 
          lmat = lmat, 
          lambda = 10^-8)
})