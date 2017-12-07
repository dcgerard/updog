context("test util_funs")

test_that("xi_fun works", {
  p <- seq(0, 1, length = 20)
  eps <- 0.1
  h <- 0.2

  eta <- p * (1 - eps) + (1 - p) * eps
  xi <- eta / (h * (1 - eta) + eta)

  expect_equal(xi_fun(p = p, eps = eps, h = h), xi)
}
)
