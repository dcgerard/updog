context("flexdog")

test_that("flexdog works", {
  refvec    <- c(1, 2, 3)
  sizevec   <- c(6, 5, 4)
  ploidy    <- 6
  model     <- c("ash", "flex")
  verbose   <- TRUE
  mean_bias <- 0
  var_bias  <- 1
  mean_seq  <- -4.7
  var_seq   <- 1
  seq       <- 0.005
  bias      <- 1
  od        <- 0.001
  mode      <- NULL
}
)


test_that("get_probk_vec works", {
  pivec <- c(0.3, 0.2, 0.5)

  expect_equal(get_probk_vec(pivec = pivec, model = "flex", mode = -1), pivec)

  expect_equal(get_probk_vec(pivec = pivec, model = "ash", mode = 2.5), c(0.1, 0.2, 0.7))

  expect_equal(get_probk_vec(pivec = pivec, model = "ash", mode = 1.5), c(0.15, 0.35, 0.5))

  expect_equal(get_probk_vec(pivec = pivec, model = "ash", mode = 0.5), c(0.3, 0.45, 0.25))

  expect_equal(get_probk_vec(pivec = pivec, model = "ash", mode = -0.5), c(0.4 + 0.5/3, 0.1 + 0.5/3, 0.5/3))


}
)
