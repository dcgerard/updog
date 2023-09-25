test_that("stan simplex works", {
  expect_equal(real_to_simplex(y = rep(0, 4)), rep(1/5, 5))

  set.seed(1)
  expect_equal(sum(real_to_simplex(y = rnorm(10))), 1)
})
