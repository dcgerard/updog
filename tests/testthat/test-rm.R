test_that("stan simplex works", {
  expect_equal(c(real_to_simplex(y = rep(0, 4))), rep(1/5, 5))

  expect_equal(c(real_to_simplex(c(2, 3, 1))), real_to_simplex_r(c(2, 3, 1)))

  set.seed(1)
  expect_equal(sum(c(real_to_simplex(y = rnorm(10)))), 1)
})

test_that("c obj function works", {
  wvec <- 1:5
  y <- c(-1, 2)
  expect_equal(
    obj_rm(y = y, weight_vec = wvec),
    obj_rm_r(weight_vec = wvec, gam = y)
  )
})

test_that("dreal_to_simplex() works", {
  y <- c(0, 1)
  myenv <- new.env()
  assign(x = "y", value = y, envir = myenv)
  nout <- stats::numericDeriv(
    quote(
      real_to_simplex(y = y)
      ),
    "y",
    myenv
    )
  expect_equal(
    attr(nout, "gradient"),
    dreal_to_simplex_dy(y = y)
  )
})

test_that("dq_dp works", {
  qfun <- function(p) stats::convolve(p, rev(p), type = "open")

  p <- c(0.1, 0.5, 0.4)
  myenv <- new.env()
  assign(x = "p", value = p, envir = myenv)
  nout <- stats::numericDeriv(
    quote(
      qfun(p = p)
    ),
    "p",
    myenv
  )

  expect_equal(
    dq_dp(p = p),
    attr(nout, "gradient"),
    tolerance = 10^-6
  )
})

test_that("dobjrm_dy works", {
  wvec <- 1:5
  y <- c(-1, 2)

  myenv <- new.env()
  assign(x = "y", value = y, envir = myenv)
  assign(x = "weight_vec", value = wvec, envir = myenv)
  nout <- stats::numericDeriv(
    quote(
      obj_rm(y = y, weight_vec = weight_vec)
    ),
    "y",
    myenv
  )

  expect_equal(
    c(dobjrm_dy(y = y, weight_vec = wvec)),
    c(attr(nout, "gradient")),
    tolerance = 10^-7
  )
})

test_that("finite rm values", {
  weight_vec <- c(4.92884695402454, 6.00089833550996, 12.5633447348355, 1.50690996815318, 7.47681185216578e-09)
  gam <- c(0, 0)

  expect_error(
    {
      optim_out <- stats::optim(
          par = gam,
          fn = obj_rm,
          gr = dobjrm_dy,
          method = "L-BFGS-B",
          weight_vec = weight_vec,
          control = list(fnscale = -1)
      )
      },
  NA)

  weight_vec <- c(1.02560280389489e-15, 2.22988144615752, 22.4813069872564, 0.288811566586071, 4.2340100687496e-16)
  gam <- c(0, 0)

  obj_rm(weight_vec = weight_vec, gam = gam)

  optim_out <- stats::optim(
    par = gam,
    fn = obj_rm,
    gr = drmlike_dq,
    method = "L-BFGS-B",
    weight_vec = weight_vec,
    control = list(fnscale = -1)
  )

  phwep <- rm_em(weight_vec = weight_vec, pvec = c(1/3, 1/3, 1/3))
  pgrad <- real_to_simplex(optim_out$par)
  expect_equal(phwep, pgrad, tol = 0.001)
  expect_equal(
    rm_llike(weight_vec = weight_vec, pvec = phwep),
    rm_llike(weight_vec = weight_vec, pvec = pgrad),
    tol = 1e-5
  )

  # bound <- 20
  # ploidy <- 3
  # real_to_simplex(y = rep(bound, ploidy - 1))

})

test_that("numeric derivative works for test case", {
  y <- c(-5.754, 46.38)
  p <- real_to_simplex(y)
  q <- stats::convolve(p, rev(p), type = "open")
  weight_vec <- c(1.02560280389489e-15, 2.22988144615752, 22.4813069872564, 0.288811566586071, 4.2340100687496e-16)

  j_dp_dy <- dreal_to_simplex_dy(y = y)
  j_dq_dp <- dq_dp(p = p)
  j_df_dq <- drmlike_dq(q = q, weight_vec = weight_vec)
  ## dobjrm_dy(y = y, weight_vec = weight_vec)

  expect_true(all(!is.nan(t(j_dp_dy) %*% t(j_dq_dp) %*% j_df_dq)))

})

