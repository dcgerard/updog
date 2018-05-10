context("util funs")

test_that("xi_fun and eta_fun works", {
  p <- seq(0, 1, length = 20)
  eps <- 0.1
  h <- 0.2
  eta <- p * (1 - eps) + (1 - p) * eps
  expect_equal(eta, eta_fun(p, eps))
  xi <- eta / (h * (1 - eta) + eta)
  expect_equal(xi_fun(p = p, eps = eps, h = h), xi)

  eps <- seq(0.1, 0.9, length = 20)
  h   <- seq(0.1, 0.9, length = 20)
  eta <- p * (1 - eps) + (1 - p) * eps
  expect_equal(eta, eta_fun(p, eps))
  xi <- eta / (h * (1 - eta) + eta)
  expect_equal(xi_fun(p = p, eps = eps, h = h), xi)

  # microbenchmark::microbenchmark(
  #   p * (1 - eps) + (1 - p) * eps,
  #   eta_fun(p, eps)
  # )
}
)

test_that("log_sum_exp works", {
  set.seed(1)
  x <- abs(rnorm(10))

  naive_way <- log(sum(exp(x)))
  expect_equal(naive_way, log_sum_exp(x))

  expect_equal(log_sum_exp(c(-Inf, 10)), 10)
  expect_equal(log_sum_exp(c(-Inf, -Inf)), -Inf)
  expect_equal(log_sum_exp(10), 10)
  expect_equal(log_sum_exp(-Inf), -Inf)
}
)

test_that("logit works", {
  expect_equal(logit(0.1), log(0.1 / 0.9))

  expect_equal(expit(logit(0.1)), 0.1)
  expect_equal(logit(expit(-8)), -8)
}
)

test_that("log_sum_exp and log_sum_exp_2 are equal", {
  expect_equal(log_sum_exp(c(0.1, 0.2)),
               log_sum_exp_2(0.1, 0.2))
}
)

test_that("log_sum_exp and log_sum_exp_2 works with -Inf", {
  expect_true(!is.nan(log_sum_exp(c(-Inf, -Inf))))
  expect_true(!is.nan(log_sum_exp_2(-Inf, -Inf)))
  expect_true(is.infinite(log_sum_exp(c(-Inf, -Inf))))
  expect_true(is.infinite(log_sum_exp_2(-Inf, -Inf)))
})
