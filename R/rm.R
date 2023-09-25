## Random mating functions

#' Objective used for random mating prior
#'
#' @noRd
obj_rm <- function(weight_vec, gam) {
  pvec <- real_to_simplex(y = gam)
  qvec <- stats::convolve(pvec, rev(pvec), type = "open")
  ## convolve can produce slightly negative or slightly above 1 values based
  ## on numerical precision.
  qvec[qvec < 0] <- 0
  qvec[qvec > 1] <- 1
  qvec <- qvec / sum(qvec)
  lvec <- weight_vec * log(qvec)

  ## remove where weight_vec is 0,
  ## because qvec = 0 and weight_vec = 0 should not contribute
  lvec[weight_vec < sqrt(.Machine$double.eps)] <- 0
  return(sum(lvec))
}

#' stan parameterization of simplex
#'
#' \url{https://mc-stan.org/docs/2_19/reference-manual/simplex-transform-section.html}
#'
#' @noRd
real_to_simplex <- function(y) {
  expit <- function(x) 0.5 + 0.5 * tanh(x / 2)

  K <- length(y)
  z <- expit(y + log(1 / (K + 1 - seq_along(y))))

  x <- rep(NA_real_, length.out = K + 1)
  x[[1]] <- z[[1]]

  if (K == 1) {
    x[[2]] <- 1 - x[[1]]
    return(x)
  }

  for (i in 2:K) {
    x[[i]] <- (1 - sum(x[1:(i-1)])) * z[[i]]
  }

  x[[K+1]] <- 1 - sum(x[1:K])

  return(x)
}

