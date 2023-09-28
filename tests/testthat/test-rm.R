test_that("stan simplex works", {
  expect_equal(real_to_simplex(y = rep(0, 4)), rep(1/5, 5))

  set.seed(1)
  expect_equal(sum(real_to_simplex(y = rnorm(10))), 1)
})

test_that("finite rm values", {
  weight_vec <- c(4.92884695402454, 6.00089833550996, 12.5633447348355, 1.50690996815318, 7.47681185216578e-09)
  gam <- c(0, 0)

  expect_error(
    {
      optim_out <- stats::optim(
          par = gam,
          fn = obj_rm,
          method = "L-BFGS-B",
          lower = c(-10, -10), upper = c(10, 10),
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
    method = "L-BFGS-B",
    lower = c(-10, -10), upper = c(10, 10),
    weight_vec = weight_vec,
    control = list(fnscale = -1)
  )

  # bound <- 20
  # ploidy <- 3
  # real_to_simplex(y = rep(bound, ploidy - 1))

})

