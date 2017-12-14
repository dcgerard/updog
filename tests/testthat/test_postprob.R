context("objective function")

test_that("post_prob works", {

  cppway <- post_prob(2, 6, -1, 1, 0.5, 0.01)
  rway <- pnorm(qnorm(pbetabinom(2, 6, 0.5, 0.01, FALSE)), mean = -1, sd = 1) -
    pnorm(qnorm(pbetabinom(1, 6, 0.5, 0.01, FALSE)), mean = -1, sd = 1)

  expect_equal(cppway, rway)
}
)

test_that("pen_bias works", {
  cppway <- pen_bias(h = 0.8, mu_h = 1, sigma2_h = 1)
  rway   <- -log(0.8) - (log(0.8) - 1) ^ 2 / (2 * 1)
  expect_equal(cppway, rway)
}
)

test_that("pen_seq_error works", {
  cppway <- pen_seq_error(0.2, -4, 1)
  rway   <- -log(0.2 * 0.8) - (log(0.2 / 0.8) + 4) ^ 2 / (2 * 1)
  expect_equal(cppway, rway)
}
)

test_that("compute_all_log_bb returns NA", {

  refmat <- matrix(c(1,2,3,NA), nrow = 2)
  sizemat <- matrix(c(NA,3,4,5), nrow = 2)
  ploidy <- 4
  seq <- c(0.2, 0.1)
  bias <- c(1, 1.3)
  od <- c(0.01, 0.02)

  bbdenseout <- compute_all_log_bb(refmat = refmat, sizemat = sizemat, ploidy = ploidy, seq = seq, bias = bias, od = od)

  expect_true(all(is.na(bbdenseout[1, 1, ])))
  expect_true(all(is.na(bbdenseout[2, 2, ])))
  expect_false(all(is.na(bbdenseout[1, 2, ])))
  expect_false(all(is.na(bbdenseout[2, 1, ])))
}
)
