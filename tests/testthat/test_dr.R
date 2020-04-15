context("Double Reduction")

test_that("get_bivalent_probs_dr works", {
  ploidy <- 8
  blist <- get_bivalent_probs(ploidy = ploidy)
  drlist <- get_bivalent_probs_dr(ploidy = ploidy)

  expect_equal(drlist$probmat[drlist$penvec, ], blist$probmat)
  expect_equal(drlist$pmat[drlist$penvec, ], blist$pmat)
  expect_equal(drlist$lvec[drlist$penvec], blist$lvec)
})
