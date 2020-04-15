context("F1 genotyping")

test_that("get intuitive F1 genotyping", {
  load("testdat.RData")
  ploidy <- 2

  fout <- flexdog_full(refvec = refvec,
                       sizevec  = sizevec,
                       ploidy   = ploidy,
                       model    = "f1",
                       p1ref    = p1ref,
                       p1size   = p1size,
                       p2ref    = p2ref,
                       p2size   = p2size)
  expect_equal(fout$par$p1geno, 0)
  expect_equal(fout$par$p2geno, 1)
})


test_that("flex_update_pivec gives intuitive results", {
  weight_vec <- c(70.87955, 64.10825, 0.01219)
  control <- list(qarray = structure(c(1, 0.5, 0, 0.5, 0.25, 0, 0, 0, 0, 0,
                                       0.5, 1, 0.5, 0.5, 0.5, 1, 0.5, 0, 0,
                                       0, 0, 0, 0.25, 0.5, 0, 0.5, 1),
                                     .Dim = c(3L, 3L, 3L),
                                     .Dimnames = list(parent1 = c("aa", "Aa", "AA"),
                                                      parent2 = c("aa", "Aa", "AA"),
                                                      offspring = c("aa", "Aa", "AA"))),
                  fs1_alpha = 0.001,
                  outliers = FALSE,
                  pivec = c(0.501248485507885, 0.413481115690214, 0.0852703988019006),
                  p1_lbb = c(-0.0880572530585763, -26.8647138792035, -78.4850500935916),
                  p2_lbb = c(-40.6462413090379, -5.94237778106282, -54.0230897301435))

  pout <- flex_update_pivec(weight_vec = weight_vec, model = "f1", control = control)
  expect_equal(pout$par$p1geno, 0)
  expect_equal(pout$par$p2geno, 1)
})
