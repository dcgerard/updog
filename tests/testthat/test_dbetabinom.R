test_that("dbetabinom works", {

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
  probseq <- c(0, 0.2, 0.25, 0.4, 0.5, 0.7, 0.7, 1)
  rway <- qbinom(p = probseq, size = 2, prob = 0.5)

  myway <- rep(NA, length = length(probseq))
  for (index in seq_along(probseq)) {
    myway[index] <- qbetabinom_double(p = probseq[index], size = 2, mu = 0.5, rho = 0)
  }

  myway2 <- qbetabinom(p = probseq, size = 2, mu = 0.5, rho = 0)

  expect_equal(rway, myway)
  expect_equal(rway, myway2)
}
)


test_that("rbetabinom works", {
  set.seed(1)
  rbbout <- c(table(rbetabinom(n = 10000, size = 2, mu = 0.3, rho = 0.1))) / 10000
  dbbout <- dbetabinom(x = 0:2, size = 2, mu = 0.3, rho = 0.1, log = FALSE)
  expect_equal(dbbout, rbbout, tolerance = 0.01, ignore_attr = TRUE)
}
)
