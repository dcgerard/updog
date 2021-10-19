context("rflexdog")

test_that("rgeno works", {
  skip_on_os(os = "mac", arch = "aarch64")

  pivec <- runif(7)
  pivec <- pivec / sum(pivec)
  trash <- rgeno(n = 10, ploidy = 6, model = "hw", allele_freq = 0.5)
  trash <- rgeno(n = 10, ploidy = 6, model = "f1", p1geno = 1, p2geno = 1)
  trash <- rgeno(n = 10, ploidy = 6, model = "s1", p1geno = 1)
  trash <- rgeno(n = 10, ploidy = 6, model = "flex", pivec = pivec)
  trash <- rgeno(n = 10, ploidy = 6, model = "uniform")
  trash <- rgeno(n = 10, ploidy = 6, model = "bb", allele_freq = 0.5, od = 0.1)
  trash <- rgeno(n = 10, ploidy = 6, model = "norm", mu = 1, sigma = 1)
  expect_error(rgeno(n = 10, ploidy = 6, model = "hw", allele_freq = 0.5, p1geno = 1))
})

test_that("rflexdog works", {
  skip_on_os(os = "mac", arch = "aarch64")

  set.seed(1)
  sizevec <- rep(1000000, length = 7)
  geno    <- 0:6
  ploidy  <- 6
  seq     <- 0
  bias    <- 1
  od      <- 0
  tvec <- rflexdog(sizevec = sizevec, geno = geno, ploidy = ploidy, seq = seq, bias = bias, od = od) / 1000000
  expect_equal(tvec, (0:6) / 6, tol = 10^-2)
})
