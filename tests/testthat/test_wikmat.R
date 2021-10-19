context("wik_mat and probvec")

test_that("wik_mat gives probvec when sizevec == 0", {
  skip_on_os(os = "mac", arch = "aarch64")

  refvec <- c(1, 1, 0)
  sizevec <- c(1, 2, 0)
  ploidy <- 4
  probk_vec <- dbinom(x = 0:ploidy, size = ploidy, prob = 0.5)
  seq <- 0.01
  bias <- 1
  od <- 0.01

  wik_mat <- get_wik_mat(probk_vec = probk_vec,
                         refvec    = refvec,
                         sizevec   = sizevec,
                         ploidy    = ploidy,
                         seq       = seq,
                         bias      = bias,
                         od        = od)

  expect_equal(wik_mat[3, ], probk_vec)
})
