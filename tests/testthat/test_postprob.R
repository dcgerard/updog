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
