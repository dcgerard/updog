## Random mating functions



#' Objective used for random mating prior on simplex scale
#'
#' @param weight_vec vector of weights
#' @param pvec vector of gamete frequencies
#'
#' @noRd
rm_llike <- function(weight_vec, pvec) {
  qvec <- stats::convolve(pvec, rev(pvec), type = "open")
  qvec[qvec < 0] <- 0
  qvec[qvec > 1] <- 1
  qvec <- qvec / sum(qvec)
  lvec <- weight_vec * log(qvec)
  lvec[weight_vec < sqrt(.Machine$double.eps)] <- 0
  return(sum(lvec))
}

#' Objective used for random mating prior on unconstrained scale
#'
#' @param weight_vec vector of weights
#' @param gam vector of unconstrained parameterization on gamete frequency simplex
#'
#' @noRd
obj_rm_r <- function(weight_vec, gam) {
  pvec <- real_to_simplex(y = gam)
  return(rm_llike(weight_vec = weight_vec, pvec = pvec))
}

#' stan parameterization of simplex
#'
#' \url{https://mc-stan.org/docs/2_19/reference-manual/simplex-transform-section.html}
#'
#' @noRd
real_to_simplex_r <- function(y) {
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

#' EM algorithm to maximize rm_llike
#'
#' Slower than gradient ascent, but probably more stable.
#'
#' @param weight_vec vector of weights
#' @param pvec initialization
#' @param tol tolerance for stopping
#' @param maxit maximum number of iterations
#'
#' @noRd
rm_em <- function(weight_vec,
                  pvec,
                  tol = 10^-3,
                  maxit = 100) {
  ploidy <- length(weight_vec) - 1
  ll <- rm_llike(weight_vec = weight_vec, pvec = pvec)

  ## Initialize parameters -----
  paramdf <- as.data.frame(which(upper.tri(matrix(nrow = ploidy / 2 + 1, ncol = ploidy / 2 + 1), diag = TRUE), arr.ind = TRUE))
  names(paramdf) <- c("i", "j")
  paramdf$geno <- paramdf$i + paramdf$j - 2
  paramdf$w <- NA_real_
  paramdf$xi <- NA_real_

  etavec <- rep(NA_real_, length = ploidy / 2 + 1)

  ## Run EM algorthm
  i <- 1
  err <- Inf
  while (i < maxit && err > tol) {
    llold <- ll
    ## One fixed point iteration ----
    paramdf$w <- pvec[paramdf$i] * pvec[paramdf$j]
    paramdf$w[paramdf$i != paramdf$j] <- 2 * paramdf$w[paramdf$i != paramdf$j]
    sumout <- by(data = paramdf, INDICES = paramdf$geno, FUN = function(x) sum(x$w), simplify = TRUE)
    paramdf$w <- paramdf$w / as.vector(sumout)[match(paramdf$geno, names(sumout))]
    paramdf$xi <- paramdf$w * weight_vec[paramdf$geno + 1]
    paramdf$xi[paramdf$i == paramdf$j] <- 2 * paramdf$xi[paramdf$i == paramdf$j]

    for (j in seq_len(ploidy / 2 + 1)) {
      etavec[[j]] <- sum(paramdf$xi[paramdf$i == j | paramdf$j == j])
    }

    ## normalize ----
    pvec <- etavec / sum(etavec)

    ## calculate log-likelihood ----
    ll <- rm_llike(weight_vec = weight_vec, pvec = pvec)

    if (ll - llold < -10^-6) {
      stop("rmem: log-likelihood is not increasing")
    }

    err <- ll - llold

    i <- i + 1
  }

  return(pvec)
}
