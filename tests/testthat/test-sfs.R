test_that("SAF calculation", {
  gl <- structure(c(-128.496694754982, -37.6997463101416, -91.1772057718566,
-75.3587315691724, -43.6295062153546, -42.6684543107057, -2.96923496083154,
-20.7156088526214, -11.8886878076436, -3.6030551107411, -12.3400919517075,
-13.0517667422853, -3.86245366217879, -3.55892819431222, -8.50466736634985,
-2.9329228265788, -42.1503475046044, -8.31309900548329, -18.6469972318237,
-31.8809740741026, -31.1618892319976, -113.684525143024, -53.7832861175473,
-78.3170711013334, -96.5512806199282), dim = c(5L, 5L), dimnames = list(
    ind = c("PotVar0000693", "PotVar0003381", "PotVar0004919",
    "PotVar0009402", "PotVar0009518"), variable = c("logL_0",
    "logL_1", "logL_2", "logL_3", "logL_4")))

  saf1 <- saf_naive(gl)
  saf2 <- saf_ll(gl)

  expect_equal(saf1, saf2)
})


test_that("matrix convolution works", {
  A <- matrix(1:12, nrow = 4)

  c1 <- convolve_m(A[1:2, ])
  c2 <- convolve(A[1, ], rev(A[2, ]), type = "open")
  expect_equal(c1, c2)

  c3 <- convolve_m(A[1:3, ])
  c4 <- convolve(c2, rev(A[3, ]), type = "open")
  expect_equal(c3, c4)

  c5 <- convolve_m(A[1:4, ])
  c6 <- convolve(c4, rev(A[4, ]), type = "open")
  expect_equal(c5, c6)

})
