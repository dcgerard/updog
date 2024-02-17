test_that("stan simplex works", {
  expect_equal(c(real_to_simplex(y = rep(0, 4))), rep(1/5, 5))

  expect_equal(c(real_to_simplex(c(2, 3, 1))), real_to_simplex_r(c(2, 3, 1)))

  set.seed(1)
  expect_equal(sum(c(real_to_simplex(y = rnorm(10)))), 1)
})

test_that("c obj function works", {
  wvec =1:5
  y =c(-1, 2)
  expect_equal(
    obj_rm(y = y, weight_vec = wvec),
    obj_rm_r(weight_vec = wvec, gam = y)
  )
})

test_that("dreal_to_simplex() works", {
  y =c(0, 1)
  myenv =new.env()
  assign(x = "y", value = y, envir = myenv)
  nout =stats::numericDeriv(
    quote(
      real_to_simplex(y = y)
      ),
    "y",
    myenv
    )
  expect_equal(
    attr(nout, "gradient"),
    dreal_to_simplex_dy(y = y)
  )
})

test_that("dq_dp works", {
  qfun =function(p) stats::convolve(p, rev(p), type = "open")

  p =c(0.1, 0.5, 0.4)
  myenv =new.env()
  assign(x = "p", value = p, envir = myenv)
  nout =stats::numericDeriv(
    quote(
      qfun(p = p)
    ),
    "p",
    myenv
  )

  expect_equal(
    dq_dp(p = p),
    attr(nout, "gradient"),
    tolerance = 10^-6
  )
})

test_that("dobjrm_dy works", {
  wvec =1:5
  y =c(-1, 2)

  myenv =new.env()
  assign(x = "y", value = y, envir = myenv)
  assign(x = "weight_vec", value = wvec, envir = myenv)
  nout =stats::numericDeriv(
    quote(
      obj_rm(y = y, weight_vec = weight_vec)
    ),
    "y",
    myenv
  )

  expect_equal(
    c(dobjrm_dy(y = y, weight_vec = wvec)),
    c(attr(nout, "gradient")),
    tolerance = 10^-7
  )
})

test_that("finite rm values", {
  weight_vec =c(4.92884695402454, 6.00089833550996, 12.5633447348355, 1.50690996815318, 7.47681185216578e-09)

  expect_no_error({
    em1 <- rm_em(weight_vec = weight_vec, pvec = rep(1/3, 3))
  })

  weight_vec <- c(0, 0, 0, 0, 10)
  pvec <- c(0, 0, 1)
  expect_true(!all(is.nan(rm_em(weight_vec = weight_vec, pvec = pvec))))

  ## Try out corner cases ----
  wmat <- as.matrix(
    expand.grid(w0 = c(0, 1),
                w1 = c(0, 1),
                w2 = c(0, 1),
                w3 = c(0, 1),
                w4 = c(0, 1))
  )

  for (i in seq_len(nrow(wmat))) {
    weight_vec <- wmat[i, ]
    expect_true(!any(is.nan(rm_em(weight_vec = weight_vec, pvec = rep(1/3, 3)))))
  }
})

test_that("EM and gradient ascient work on easy case", {
  weight_vec <- c(4, 5, 2, 6, 7)
  gam <- c(0, 0)
  pvec <- rep(1/3, 3)


  emout <- rm_em(weight_vec = weight_vec, pvec = pvec)

  gradout <- stats::optim(
    par = gam,
    fn = obj_rm,
    gr = dobjrm_dy,
    method = "L-BFGS-B",
    control = list(fnscale = -1),
    weight_vec = weight_vec
  )


  expect_equal(
    real_to_simplex(gradout$par),
    emout,
    tolerance = 1e-3
  )
})

test_that("Difficult SNPs work", {
  load("./rmdat.RData")
  expect_no_error(
    mout <- multidog(refmat = refmat, sizemat = sizemat, ploidy = 4, model = "rm")
  )
})

test_that("rm prior is symmetric", {
  ploidy <- 4L
  genovec1 <- c(2L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 2L, 1L, 2L)
  genovec2 <- 4L - genovec1
  refvec1 <- c(21L, 0L, 0L, 1L, 0L, 1L, 0L, 14L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 12L, 0L, 0L, 1L, 22L, 12L, 16L)
  refvec2 <- 45L - refvec1
  sizevec <- rep(45L, 25)

  uout1 <- flexdog(refvec = refvec1, sizevec = sizevec, ploidy = ploidy, model = "rm", var_bias = Inf)
  uout2 <- flexdog(refvec = refvec2, sizevec = sizevec, ploidy = ploidy, model = "rm", var_bias = Inf)

  4 - uout2$geno
  uout1$geno

  log(uout1$bias)
  log(uout2$bias)
})
