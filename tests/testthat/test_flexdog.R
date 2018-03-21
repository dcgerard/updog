context("flexdog")

test_that("flexdog works", {
  refvec    <- 1:20
  sizevec   <- 40:21

  refvec[2]  <- NA
  sizevec[3] <- NA

  # snpdat <- updog::snpdat
  # library(tidyverse)
  # snpdat %>% filter(snp == "SNP1") ->
  #   smalldat
  # refvec <- smalldat$counts[!is.na(smalldat$counts)]
  # sizevec <- smalldat$size[!is.na(smalldat$size)]

  ploidy    <- 4
  model     <- "ash"
  verbose   <- TRUE
  mean_bias <- 0
  var_bias  <- 1
  mean_seq  <- -4.7
  var_seq   <- 1
  seq       <- 0.005
  bias      <- 1
  od        <- 0.001
  mode      <- NULL
  itermax   = 10
  tol       = 10^-2

  fout <- flexdog(refvec = refvec, sizevec = sizevec,
                  ploidy = ploidy, model = "hw")
  fout <- flexdog(refvec = refvec, sizevec = sizevec,
                  ploidy = ploidy, model = "flex")
  pl <- plot(fout)
}
)

test_that("get_uni_rep is inverse of get_probk_vec", {
  probvec <- dbinom(0:6, 6, 0.4)
  pout <- get_uni_rep(probvec)
  probvec2 <- get_probk_vec(pivec = pout$pivec,
                            model = "ash",
                            mode = pout$mode - 1)
  expect_equal(probvec, probvec2)

  probvec <- dbinom(0:6, 6, 0.9)
  pout <- get_uni_rep(probvec)
  probvec2 <- get_probk_vec(pivec = pout$pivec,
                            model = "ash",
                            mode = pout$mode - 1)
  expect_equal(probvec, probvec2)
})

test_that("get right mode in initialization", {
  ploidyvec <- 2:30
  for (ploidy in ploidyvec) {
    modevec <- 0:(ploidy - 1) + 0.5
    for (index in 1:length(modevec)) {
      mode <- modevec[index]
      pvec_init <- stats::dbinom(x = 0:ploidy, size = ploidy,
                                 prob = floor(mode) / ploidy)
      pout <- get_uni_rep(pvec_init)
      expect_equal(pout$mode - 1, mode)
    }
  }
})

test_that("get_probk_vec works", {
  pivec <- c(0.3, 0.2, 0.5)

  expect_equal(get_probk_vec(pivec = pivec, model = "flex", mode = -1), pivec)

  expect_equal(get_probk_vec(pivec = pivec, model = "ash", mode = 2.5), c(0.1, 0.2, 0.7))

  expect_equal(get_probk_vec(pivec = pivec, model = "ash", mode = 1.5), c(0.15, 0.35, 0.5))

  expect_equal(get_probk_vec(pivec = pivec, model = "ash", mode = 0.5), c(0.3, 0.45, 0.25))

  expect_equal(get_probk_vec(pivec = pivec, model = "ash", mode = -0.5), c(0.4 + 0.5/3, 0.1 + 0.5/3, 0.5/3))


}
)
