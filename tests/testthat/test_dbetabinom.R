context("dbetabinom")

test_that("dbetabinom works", {

  if (requireNamespace("SuppDists", quietly = TRUE)) {

    dbetabinom_suppdists <- function(x, size, mu, rho, log) {
      alpha <- mu * (1 - rho) / rho
      beta  <- (1 - mu) * (1 - rho) / rho
      SuppDists::dghyper(x = x, a = -alpha, k = size, N = -alpha - beta, log = log)
    }

    expect_equal(dbetabinom_suppdists(1, 2, 0.5, 0.2, TRUE),
                 dbetabinom(1, 2, 0.5, 0.2, TRUE))

    # microbenchmark::microbenchmark(
    #   dbetabinom_suppdists(1, 2, 0.5, 0.2, TRUE),
    #   mupdog::dbetabinom(1, 2, 0.5, 0.2, TRUE),
    #   VGAM::dbetabinom(1, 2, 0.5, 0.2, TRUE)
    # )
  }

  trash <- dbetabinom(c(1, 2), c(2, 3), 0.5, 0.2, TRUE)
  trash <- dbetabinom(c(1, 2), c(2, 3), c(0.5, 0.2), 0.2, TRUE)
  trash <- dbetabinom(c(1, 2), c(2, 3), 0.5, c(0.5, 0.2), TRUE)
  trash <- dbetabinom(c(1, 2), c(2, 3), c(0.5, 0.2), c(0.5, 0.2), TRUE)

  trash <- capture.output(expect_error(dbetabinom(3, 2, 0.5, 0.2, TRUE)))
  trash <- capture.output(expect_error(dbetabinom(-1, 2, 0.5, 0.2, TRUE)))

  expect_equal(dbetabinom(0, 2, 0, 0.5, FALSE), 1)
  expect_equal(dbetabinom(2, 2, 1, 0.5, FALSE), 1)
  expect_equal(dbetabinom(2, 2, 0, 0.5, FALSE), 0)
  expect_equal(dbetabinom(0, 2, 1, 0.5, FALSE), 0)

  expect_equal(dbetabinom(0, 2, 0.2, 1, FALSE), 0.8)
  expect_equal(dbetabinom(1, 2, 0.2, 1, FALSE), 0)
  expect_equal(dbetabinom(2, 2, 0.2, 1, FALSE), 0.2)

  expect_equal(dbinom(x = 1, size = 2, prob = 0.2, log = FALSE),
               dbetabinom(1, 2, 0.2, 0, FALSE))

  trash <- capture.output(expect_error(dbetabinom(1, 2, -1, 0.3, TRUE)))
  trash <- capture.output(expect_error(dbetabinom(1, 2, 0.5, -1, TRUE)))

}
)

test_that("pbetabinom works ok", {

  if (requireNamespace("SuppDists", quietly = TRUE)) {

    pbetabinom_suppdists <- function(q, size, mu, rho, log) {
      alpha <- mu * (1 - rho) / rho
      beta  <- (1 - mu) * (1 - rho) / rho
      SuppDists::pghyper(q = q, a = -alpha, k = size, N = -alpha - beta, log = log)
    }

    expect_equal(pbetabinom_suppdists(4, 6, 0.5, 0.2, TRUE),
                 pbetabinom(4, 6, 0.5, 0.2, TRUE))

    # microbenchmark::microbenchmark(
    #   pbetabinom_suppdists(4, 6, 0.5, 0.2, TRUE),
    #   pbetabinom(4, 6, 0.5, 0.2, TRUE),
    #   VGAM::pbetabinom(4, 6, 0.5, 0.2, TRUE)
    #   )
  }

  trash <- pbetabinom(c(3, 4), c(6, 6), 0.5, 0.2, TRUE)
  trash <- pbetabinom(c(3, 4), c(6, 6), c(0.5, 0.7), 0.2, TRUE)
  trash <- pbetabinom(c(3, 4), c(6, 6), c(0.5, 0.7), c(0.5, 0.2), TRUE)
  trash <- pbetabinom(c(3, 4), c(6, 6), 0.5, c(0.5, 0.2), TRUE)

  expect_equal(pbetabinom(-1, 2, 0.5, 0.5, FALSE), 0)
  expect_equal(pbetabinom(3, 2, 0.5, 0.5, FALSE), 1)

  expect_equal(pbetabinom(0, 2, 0, 0.5, FALSE), 1)
  expect_equal(pbetabinom(1, 2, 1, 0.5, FALSE), 0)

}
)
