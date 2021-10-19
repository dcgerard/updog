context("flexdog")

test_that("flexdog works", {
  skip_on_os(os = "mac", arch = "aarch64")

  refvec    <- 1:20
  sizevec   <- 40:21

  refvec[2]  <- NA
  sizevec[3] <- NA

  # data("snpdat")
  # library(tidyverse)
  # snpdat %>% filter(snp == "SNP1") ->
  #   smalldat
  # refvec <- smalldat$counts[!is.na(smalldat$counts)]
  # sizevec <- smalldat$size[!is.na(smalldat$size)]

  ploidy      <- 4
  model       <- "bb"
  verbose     <- TRUE
  mean_bias   <- 0
  var_bias    <- 1
  mean_seq    <- -4.7
  var_seq     <- 1
  seq         <- 0.005
  bias        <- 1
  od          <- 0.001
  mode        <- NULL
  itermax     <- 10
  tol         <- 10^-2
  use_cvxr    <- FALSE
  update_bias <- TRUE
  update_seq  <- TRUE
  update_od   <- TRUE
  ashpen      <- 0
  fs1_alpha   <- 10 ^ -3
  p1ref       <- NULL
  p1size      <- NULL
  p2ref       <- NULL
  p2size      <- NULL
  outliers    <- FALSE

  fout <- flexdog(refvec = refvec, sizevec = sizevec,
                  ploidy = ploidy, model = "f1", verbose = FALSE,
                  p1ref = 1, p1size = 2,
                  p2ref = 5, p2size = 10,
                  snpname = "abcdefg")
  expect_equal(fout$input$snpname, "abcdefg")
  pl <- plot(fout)
  fout <- flexdog(refvec = refvec, sizevec = sizevec,
                  ploidy = ploidy, model = "flex", verbose = FALSE)
  pl <- plot(fout)
  fout <- flexdog(refvec = refvec, sizevec = sizevec,
                  ploidy = ploidy, model = "bb", verbose = FALSE)

  fout <- flexdog(refvec = refvec, sizevec = sizevec,
                  ploidy = ploidy, model = "norm", verbose = FALSE)
  pivec <- dnorm(x = 0:ploidy, mean = fout$par$mu, sd = fout$par$sigma, log = FALSE)
  pivec <- pivec / sum(pivec)
  expect_equal(pivec, fout$gene_dist)

  expect_warning(fout <- flexdog(refvec = refvec, sizevec = sizevec,
                                 ploidy = ploidy, model = "uniform", verbose = FALSE))
  expect_equal(fout$gene_dist, rep(1 / (ploidy + 1), ploidy + 1))
  pl <- plot(fout)
  # suppressWarnings(
  # fout <- flexdog(refvec = refvec, sizevec = sizevec,
  #                 p2ref = 10, p2size = 20,
  #                 ploidy = ploidy, model = "f1ppdr", verbose = FALSE)
  # )

  flexdog(refvec = refvec,
          sizevec = sizevec,
          ploidy = ploidy,
          model = "custom",
          bias = 1,
          prior_vec = 0:ploidy / sum(0:ploidy),
          verbose = FALSE) -> fout
}
)

test_that("don't update bias, seq, od when supposed not to", {
  skip_on_os(os = "mac", arch = "aarch64")

  refvec  <- 1:20
  sizevec <- 40:21
  ploidy  <- 4
  fout <- flexdog(refvec = refvec, sizevec = sizevec,
                  ploidy = ploidy, model = "hw",
                  update_bias = FALSE, bias = 0.5,
                  verbose = FALSE)
  expect_equal(fout$bias, 0.5)

  fout <- flexdog(refvec = refvec, sizevec = sizevec,
                  ploidy = ploidy, model = "hw",
                  update_seq = FALSE, seq = 0.01,
                  verbose = FALSE)
  expect_equal(fout$seq, 0.01)

  fout <- flexdog(refvec = refvec, sizevec = sizevec,
                  ploidy = ploidy, model = "hw",
                  update_od = FALSE, od = 0.01,
                  verbose = FALSE)
  expect_equal(fout$od, 0.01)
})

test_that("fs1_alpha works", {
  skip_on_os(os = "mac", arch = "aarch64")

  refvec  <- 1:20
  sizevec <- 40:21
  ploidy  <- 4
  fout1 <- flexdog(refvec = refvec,
                   sizevec = sizevec,
                   ploidy = ploidy,
                   model = "s1",
                   fs1_alpha = 10^-4,
                   verbose = FALSE)
  expect_error(
    fout2 <- flexdog(refvec = refvec,
                     sizevec = sizevec,
                     ploidy = ploidy,
                     model = "s1",
                     fs1_alpha = "picard",
                     verbose = FALSE)
  )


  expect_equal(fout1$par$alpha, 10^-4)
})

test_that("actually using parent data", {
  skip_on_os(os = "mac", arch = "aarch64")

  refvec  <- 1:20
  sizevec <- 40:21
  ploidy  <- 6

  mcount <- 100000
  fout1 <- flexdog(refvec = refvec,
                   sizevec = sizevec,
                   ploidy = ploidy,
                   model = "s1",
                   fs1_alpha = 10^-4,
                   verbose = FALSE,
                   p1ref = mcount / ploidy,
                   p1size = mcount,
                   update_bias = FALSE,
                   bias_init = 1,
                   update_od = FALSE,
                   update_seq = FALSE)

  fout2 <- flexdog(refvec = refvec,
                   sizevec = sizevec,
                   ploidy = ploidy,
                   model = "s1",
                   fs1_alpha = 10^-4,
                   verbose = FALSE,
                   p1ref = 2 * mcount / ploidy,
                   p1size = mcount,
                   update_bias = FALSE,
                   bias_init = 1,
                   update_od = FALSE,
                   update_seq = FALSE)

  fout3 <- flexdog(refvec = refvec,
                   sizevec = sizevec,
                   ploidy = ploidy,
                   model = "s1",
                   fs1_alpha = 10^-4,
                   verbose = FALSE,
                   p1ref = 3 * mcount / ploidy,
                   p1size = mcount,
                   update_bias = FALSE,
                   bias_init = 1,
                   update_od = FALSE,
                   update_seq = FALSE)

  expect_equal(fout1$par$pgeno, 1)
  expect_equal(fout2$par$pgeno, 2)
  expect_equal(fout3$par$pgeno, 3)
})

test_that("genotype likelihoods and posteriors are consistent", {
  skip_on_os(os = "mac", arch = "aarch64")

  refvec <- c(1, 2, 1, 3, 2)
  sizevec <- c(4, 2, 2, 3, 4)
  ploidy <- 4

  fout <- flexdog(refvec = refvec,
                  sizevec = sizevec,
                  ploidy = ploidy,
                  verbose = FALSE)

  ind <- 2
  probhand <- fout$gene_dist * exp(fout$genologlike[ind,])
  probhand <- probhand / sum(probhand)
  expect_equal(
    fout$postmat[ind,],
    probhand,
    tolerance = 0.005
  )
})


