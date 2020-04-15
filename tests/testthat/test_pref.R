context("pref")

test_that("dist_from_p works", {
  expect_equal(dist_from_p(c(1,3,0), 8), c(1, 3, 3, 1, 0) / 8)
  expect_equal(dist_from_p(c(2, 1, 1), 8), c(0, 1, 1, 0, 0) / 2)

  expect_equal(dist_from_p(c(1,2,1), 8), c(0, 1, 2, 1, 0) / 4)
  expect_equal(dist_from_p(c(0,4,0), 8), c(1, 4, 6, 4, 1) / 16)
})

test_that("count_pairings works", {
  expect_equal(count_pairings(4), 3)
  expect_equal(count_pairings(6), 15)
})

test_that("count_doubles works", {
  expect_equal(count_doubles(4, 0), 1)
  expect_equal(count_doubles(4, 1), 6)
  expect_equal(count_doubles(4, 2), 3)
  expect_equal(count_doubles(4, 3), 0)


  expect_equal(count_doubles(6, 0), 1)
  expect_equal(count_doubles(6, 1), 15)
  expect_equal(count_doubles(6, 2), 45)
  expect_equal(count_doubles(6, 3), 15)
  expect_equal(count_doubles(6, 4), 0)
})

test_that("get_hyper_weights works", {
  skip("Doesn't work on windows but only during devtools::check().")
  blist <- get_bivalent_probs(6)
  hlist <- get_hyper_weights(6, 4)

  logvec <- blist$pmat[, 1] == hlist$pmat[, 1] &
    blist$pmat[, 2] == hlist$pmat[, 2] &
    blist$pmat[, 3] == hlist$pmat[, 3]

  mixdist <- colSums(blist$probmat[logvec, , drop = FALSE] * hlist$weightvec)
  hdist <- stats::dhyper(0:3, m = 4, n = 2, k = 3)

  expect_true(all(abs(mixdist - hdist) < 10 ^ -3))
})

test_that("convolve_up() works", {
  x <- c(1/4, 3/4)
  y <- c(1/8, 7/8)
  cout <- convolve_up(x, y)
  expect_equal(cout[1], x[1] * y[1])
  expect_equal(cout[2], x[2] * y[1] + x[1] * y[2])
  expect_equal(cout[3], x[2] * y[2])
  expect_equal(sum(cout), 1)
})


test_that("pp_brent_obj works", {
  firstmixweight <- 0.7
  alpha <- 0.001
  probmat <- matrix(c(0.1, 0.7, 0.2, 0.3, 0.3, 0.4), nrow = 2, byrow = TRUE)
  pvec <- c(1/6, 2/6, 3/6)
  weight_vec <- 1:5
  pp_brent_obj(firstmixweight = firstmixweight,
               probmat        = probmat,
               pvec           = pvec,
               weight_vec     = weight_vec,
               alpha          = alpha)

  f1_obj(alpha = alpha, pvec = c(0.1, 0.2, 0.3, 0.2, 0.2), weight_vec = weight_vec)
})

test_that("rgeno works with pp", {
  temp <- rgeno(n = 1,
                ploidy = 6,
                model = "f1pp",
                p1geno = 1,
                p2geno = 3,
                p1_pair_weights = 1,
                p2_pair_weights = c(0.3, 0.7))
})
