context("pref")

test_that("dist_from_p works", {
  expect_equal(dist_from_p(c(1,3,0), 8), c(1, 3, 3, 1) / 8)
  expect_equal(dist_from_p(c(2, 1, 1), 8), c(0, 1, 1, 0) / 2)

  expect_equal(dist_from_p(c(1,2,1), 8), c(0, 1, 2, 1) / 4)
  expect_equal(dist_from_p(c(0,4,0), 8), c(1, 4, 6, 4, 1) / 16)

  dist_from_p(c(0, 0, 2), 4)
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
  blist <- get_bivalent_probs(6)
  hlist <- get_hyper_weights(6, 4)

  logvec <- blist$pmat[, 1] == hlist$pmat[, 1] &
    blist$pmat[, 2] == hlist$pmat[, 2] &
    blist$pmat[, 3] == hlist$pmat[, 3]

  mixdist <- colSums(blist$probmat[logvec, , drop = FALSE] * hlist$weightvec)
  hdist <- stats::dhyper(0:3, m = 4, n = 2, k = 3)

  expect_true(all(abs(mixdist - hdist) < 10 ^ -12))
})
