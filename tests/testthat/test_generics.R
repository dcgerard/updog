context("generics")

test_that("summary.mupdog works", {
  msum <- summary.mupdog(mupout)
  byhand <- c(table(c(mupout$map_dosage[, 1], 0:mupout$input$ploidy)) - 1)
  expect_equal(msum$freq[, 1], byhand)


  byhand <- sum(mupout$postprob[1, 1, ] * 0:mupout$input$ploidy)
  expect_equal(msum$mean_dosage[1, 1], byhand)
}
)

test_that("plot.mupdog works", {
  data("mupout")
  plot.mupdog(mupout, 1)
}
)
