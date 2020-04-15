context("flexdog")

test_that("flexdog works", {
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

  fout <- flexdog(refvec = refvec, sizevec = sizevec,
                  p1ref = 10, p1size = 20,
                  ploidy = ploidy, model = "s1pp", verbose = FALSE)

  flexdog(refvec = refvec,
          sizevec = sizevec,
          ploidy = ploidy,
          model = "custom",
          bias = 1,
          prior_vec = 0:ploidy / sum(0:ploidy),
          verbose = FALSE) -> fout
}
)

test_that("f1ppdr works on SNP3", {
  skip("f1ppdr not quite yet")
  data("snpdat")
  fout <- flexdog_full(refvec  = snpdat$counts[snpdat$snp == "SNP3"],
                       sizevec = snpdat$size[snpdat$snp == "SNP3"],
                       ploidy  = 6,
                       model   = "f1ppdr",
                       verbose = FALSE, bias = 2)

  fout <- flexdog(refvec  = snpdat$counts[snpdat$snp == "SNP1"],
                       sizevec = snpdat$size[snpdat$snp == "SNP1"],
                       ploidy  = 6,
                       model   = "f1ppdr",
                       verbose = FALSE)

  fout <- flexdog(refvec  = snpdat$counts[snpdat$snp == "SNP3"],
                  sizevec = snpdat$size[snpdat$snp == "SNP3"],
                  ploidy  = 6,
                  model   = "s1pp",
                  verbose = FALSE)
})

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

test_that("don't update bias, seq, od when supposed not to", {
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


test_that("get_wik_mat_out has expected behavior", {
  ploidy    <- 4
  probk_vec <- seq(1 / (ploidy + 1), ploidy + 1)
  refvec    <- 1:5
  sizevec   <- 6:10
  seq       <- 0.001
  bias      <- 0.25
  od        <- 0.02
  wikmat <- get_wik_mat_out(probk_vec = probk_vec,
                            out_prop  = 1,
                            refvec    = refvec,
                            sizevec   = sizevec,
                            ploidy    = ploidy,
                            seq       = seq,
                            bias      = bias,
                            od        = od)
  expect_true(all(wikmat[1:5, 1:5] == 0))
  expect_true(all(wikmat[, 6] == 1))

  wikmat1 <- get_wik_mat_out(probk_vec = probk_vec,
                             out_prop  = 0,
                             refvec    = refvec,
                             sizevec   = sizevec,
                             ploidy    = ploidy,
                             seq       = seq,
                             bias      = bias,
                             od        = od)

  wikmat2 <- get_wik_mat(probk_vec = probk_vec,
                         refvec    = refvec,
                         sizevec   = sizevec,
                         ploidy    = ploidy,
                         seq       = seq,
                         bias      = bias,
                         od        = od)

  expect_equal(wikmat1[1:5, 1:5], wikmat2)
})
