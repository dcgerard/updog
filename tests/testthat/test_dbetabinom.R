context("dbetabinom")

test_that("dbetabinom works", {
  skip_on_os(os = "mac", arch = "aarch64")

  if (requireNamespace("SuppDists", quietly = TRUE)) {

    dbetabinom_suppdists <- function(x, size, mu, rho, log) {
      alpha <- mu * (1 - rho) / rho
      beta  <- (1 - mu) * (1 - rho) / rho
      SuppDists::dghyper(x = x, a = -alpha, k = size, N = -alpha - beta, log = log)
    }

    expect_equal(dbetabinom_suppdists(c(1, 2), 2, 0.5, 0.2, TRUE),
                 dbetabinom(c(1, 2), 2, 0.5, 0.2, TRUE))

    # microbenchmark::microbenchmark(
    #   dbetabinom_suppdists(1, 2, 0.5, 0.2, TRUE),
    #   updog::dbetabinom(1, 2, 0.5, 0.2, TRUE),
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
  skip_on_os(os = "mac", arch = "aarch64")

  if (requireNamespace("SuppDists", quietly = TRUE)) {

    pbetabinom_suppdists <- function(q, size, mu, rho, log) {
      alpha <- mu * (1 - rho) / rho
      beta  <- (1 - mu) * (1 - rho) / rho
      SuppDists::pghyper(q = q, a = -alpha, k = size, N = -alpha - beta, log = log)
    }

    expect_equal(pbetabinom_suppdists(c(3, 4), 6, 0.5, 0.2, TRUE),
                 pbetabinom(c(3, 4), 6, 0.5, 0.2, TRUE))

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


test_that("qbetabinom works", {
  skip_on_os(os = "mac", arch = "aarch64")

  probseq <- c(0, 0.2, 0.25, 0.4, 0.5, 0.7, 0.7, 1)
  rway <- qbinom(p = probseq, size = 2, prob = 0.5)

  myway <- rep(NA, length = length(probseq))
  for (index in seq_along(probseq)) {
    myway[index] <- qbetabinom_double(p = probseq[index], size = 2, mu = 0.5, rho = 0)
  }

  myway2 <- qbetabinom(p = probseq, size = 2, mu = 0.5, rho = 0)

  expect_equal(rway, myway, myway2)

  if (requireNamespace("SuppDists", quietly = TRUE)) {

    qbetabinom_suppdists <- function(p, size, mu, rho, log) {
      alpha <- mu * (1 - rho) / rho
      beta  <- (1 - mu) * (1 - rho) / rho
      SuppDists::qghyper(p = p, a = -alpha, k = size, N = -alpha - beta, log = log)
    }

    p0 <- dbetabinom(x = 0, size = 2, mu = 0.5, rho = 0.1, log = FALSE)
    pvec <- c(0, p0, p0 + 10^-6, 1 - 10^-6, 1)
    suppway <- qbetabinom_suppdists(pvec, 2, 0.5, 0.1, FALSE)
    myway <- rep(NA, length = length(pvec))
    for (index in seq_along(pvec)) {
      myway[index] <- qbetabinom_double(pvec[index], 2, 0.5, 0.1)
    }

    myway2 <- qbetabinom(p = pvec, size = 2, mu = 0.5, rho = 0.1)
    expect_equal(suppway, myway, myway2)


    # microbenchmark::microbenchmark(
    #   qbetabinom_suppdists(0.3, 6, 0.5, 0.2, FALSE),
    #   qbetabinom(0.3, 6, 0.5, 0.2)
    #   )
  }
}
)


test_that("rbetabinom works", {
  skip_on_os(os = "mac", arch = "aarch64")

  set.seed(1)
  rbbout <- c(table(rbetabinom(n = 10000, size = 2, mu = 0.3, rho = 0.1))) / 10000
  dbbout <- dbetabinom(x = 0:2, size = 2, mu = 0.3, rho = 0.1, log = FALSE)
  expect_equal(dbbout, rbbout, tol = 0.01, check.attributes = FALSE)
}
)
