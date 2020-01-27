context("CVXR")

# test_that("uni_em is same as CVXR on problem data.", {
#   if (requireNamespace("CVXR", quietly = TRUE)) {
#     set.seed(4)
#     cvx_prob <- readRDS(file = "cvx_prob.RDS")
#     ploidy <- cvx_prob$ploidy
#     weight_vec <- cvx_prob$weight_vec
#     lmat <- cvx_prob$lmat
#
#     ## EM way
#     pi_init <- runif(ploidy + 1)
#     pi_init <- pi_init / sum(pi_init)
#     uout <- uni_em(weight_vec = weight_vec,
#                    lmat = lmat,
#                    pi_init = pi_init,
#                    itermax = 200,
#                    obj_tol = 10^-5,
#                    lambda = .Machine$double.eps * 10)
#
#     ## CVXR way
#     cv_pi <- CVXR::Variable(1, ploidy + 1)
#     obj   <- sum(t(weight_vec) * log(cv_pi %*% lmat))
#     prob  <- CVXR::Problem(CVXR::Maximize(obj),
#                            constraints = list(sum(cv_pi) == 1,
#                                               cv_pi >= 0))
#     result <- solve(prob)
#     result$value
#     cout <- c(result$getValue(cv_pi))
#
#     expect_equal(cout, c(uout), tol = 10^-2)
#   }
# }
# )


test_that("uni_obj works when have zero weights", {
  set.seed(1)
  cvx_prob <- readRDS(file = "cvx_prob.RDS")
  weight_vec <- cvx_prob$weight_vec
  lmat <- cvx_prob$lmat
  pivec <- c(0, rep(1, length = length(weight_vec) - 1))
  pivec <- pivec / sum(pivec)
  expect_true(!is.nan(uni_obj(pivec = pivec, weight_vec = weight_vec, lmat = lmat, lambda = 0)))

  weight_vec[1] <- 0.1
  expect_true(is.infinite(uni_obj(pivec = pivec, weight_vec = weight_vec, lmat = lmat, lambda = 0)))

  ## Compare to manual calculation
  pivec <- runif(length(weight_vec))
  pivec <- pivec / sum(pivec)
  uout <- uni_obj(pivec = pivec, weight_vec = weight_vec, lmat = lmat, lambda = 0)
  mout <- sum(weight_vec * log(pivec %*% lmat))
  expect_equal(mout, uout)
}
)
