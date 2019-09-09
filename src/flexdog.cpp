#include "mupdog.h"

// functions for fitting flexdog


//' Obtain the genotype distribution given the distribution of discrete uniforms.
//'
//' @inheritParams flexdog_full
//' @param pivec The mixing probability of the i'th discrete uniform distribution.
//'
//' @author David Gerard
//'
//' @return A vector of numerics. Element k is the probability of genotype k.
//'
//' @seealso \code{\link{flexdog}} where this is used.
//'
// [[Rcpp::export]]
NumericVector get_probk_vec(NumericVector pivec, std::string model, double mode) {
  int K = pivec.length() - 1;
  NumericVector probk_vec(K + 1);
  if ((model == "flex") | (model == "hw") | (model == "f1") | (model == "s1") | (model == "uniform") |
    (model == "bb") | (model == "norm") | (model == "f1pp") | (model == "s1pp") |
    (model == "f1ppdr") | (model == "s1ppdr") | (model == "custom")) {
    probk_vec = pivec;
  } else if (model == "ash") {
    double denom; // what you divide the pi's by.
    for (int i = 0; i <= K; i++) { // iterate through pivec
      if (std::fabs((double)i - mode) < TOL) {
        probk_vec(i) += pivec(i);
      } else if (i < mode) {
        denom = (double)(std::floor(mode) - i + 1);
        for (int j = i; j <= mode; j++) { //iterate through probk_vec
          probk_vec(j) += pivec(i) / denom;
        }
      } else if (i > mode) { //iterate through probk_vec
        denom = (double)(i - std::ceil(mode) + 1);
        for (int j = i; j >= mode; j--) {
          probk_vec(j) += pivec(i) / denom;
        }
      } else {
        Rcpp::stop("get_probk_vec: How did you get here??");
      }
    }
  } else {
    Rcpp::stop("get_probk_vec: model must be one of 'flex', 'ash', 'hw', 'bb', 'norm', 'f1', 's1', 'f1pp', 's1pp', 'f1ppdr', 's1ppdr', 'uniform', or 'custom'");
  }
  return probk_vec;
}

//' Compute inner weights for updating the mixing proportions when using ash model.
//'
//' The (i,k)th element is \eqn{1(k \in F(a, i)) / |F(a,i)|}.
//'
//' @inheritParams flexdog_full
//'
//' @return A matrix of numerics. The weights used for the
//'    weighted EM algorithm in \code{\link{flexdog_full}}.
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
NumericMatrix get_inner_weights(int ploidy, double mode) {
  NumericMatrix inner_weights(ploidy + 1, ploidy + 1);
  double denom;
  for (int i = 0; i <= ploidy; i++) {
    if (std::fabs((double)i - mode) < TOL) {
      inner_weights(i, i) = 1.0;
    } else if (i < mode) {
      denom = (double)(std::floor(mode) - i + 1);
      for (int j = i; j <= mode; j++) {
        inner_weights(i, j) = 1.0 / denom;
      }
    } else if (i > mode) {
      denom = (double)(i - std::ceil(mode) + 1);
      for (int j = i; j >= mode; j--) {
        inner_weights(i, j) = 1.0 / denom;
      }
    } else {
      Rcpp::stop("get_inner_weights: How did you get here??");
    }
  }
  return(inner_weights);
}


//' E-step in \code{\link{flexdog}}.
//'
//' @inheritParams flexdog_full
//' @param probk_vec The vector of current prior probabilities of each genotype.
//'
//'
//' @return A matrix of numerics. The rows index the individuals and the
//'     columns index the genotype. These weights are used in the EM algorithm
//'     (and is indeed the E-step) in \code{\link{flexdog_full}}.
//'
//' @author David Gerard
//'
//' @seealso \code{\link{flexdog}} for the full EM algorithm.
//'
// [[Rcpp::export]]
NumericMatrix get_wik_mat(NumericVector probk_vec,
                          NumericVector refvec,
                          NumericVector sizevec,
                          int ploidy,
                          double seq,
                          double bias,
                          double od) {
  // Check input -----------------------------------------------------------
  int nind = refvec.length();
  if (nind != sizevec.length()) {
    Rcpp::stop("get_wik_mat: sizevec and refvec must have the same length.");
  }
  if (probk_vec.length() != ploidy + 1) {
    Rcpp::stop("get_wik_mat: probk_vec must have length ploidy + 1.");
  }

  // Calculate the posterior probability of each genotype -------------------
  NumericMatrix wik_mat(nind, ploidy + 1);
  NumericVector lprobk_vec = Rcpp::log(probk_vec);
  NumericVector xi(ploidy + 1);
  for (int k = 0; k <= ploidy; k++) {
    xi(k) = xi_double((double)k / (double)ploidy, seq, bias);
  }

  double sumi; // denominator to get wik for each i.
  NumericVector wvec(ploidy + 1);
  for (int i = 0; i < nind; i++) {
    for (int k = 0; k <= ploidy; k++) {
      wvec(k) = lprobk_vec(k) + dbetabinom_double(refvec(i), sizevec(i), xi(k), od, true);
    }
    sumi = log_sum_exp(wvec);
    wvec = Rcpp::exp(wvec - sumi);
    wik_mat(i, _) = wvec;
  }

  return wik_mat;
}

//' Log-likelihood that \code{\link{flexdog}} maximizes.
//'
//' @inheritParams flexdog_full
//' @param probk_vec The kth element is the prior probability of genotype k (when starting to count from 0).
//'
//' @author David Gerard
//'
//' @return The objective (marginal log-likelihood) used in
//'     \code{\link{flexdog_full}}.
//'
// [[Rcpp::export]]
double flexdog_obj(NumericVector probk_vec,
                   NumericVector refvec,
                   NumericVector sizevec,
                   int ploidy,
                   double seq,
                   double bias,
                   double od,
                   double mean_bias,
                   double var_bias,
                   double mean_seq,
                   double var_seq,
                   double mean_od,
                   double var_od) {
  // Check input -----------------------------------------------------------
  int nind = refvec.length();
  if (nind != sizevec.length()) {
    Rcpp::stop("get_wik_mat: sizevec and refvec must have the same length.");
  }
  if (probk_vec.length() != ploidy + 1) {
    Rcpp::stop("get_wik_mat: probk_vec must have length ploidy + 1.");
  }

  // Calculate the posterior probability of each genotype -------------------
  NumericVector lprobk_vec = Rcpp::log(probk_vec);
  NumericVector xi(ploidy + 1);
  for (int k = 0; k <= ploidy; k++) {
    xi(k) = xi_double((double)k / (double)ploidy, seq, bias);
  }

  // Calculate likelihood ---------------------------------------------------
  double obj = 0.0;
  NumericVector wvec(ploidy + 1);
  for (int i = 0; i < nind; i++) {
    for (int k = 0; k <= ploidy; k++) {
      wvec(k) = lprobk_vec(k) + dbetabinom_double(refvec(i), sizevec(i), xi(k), od, true);
    }
    obj = obj + log_sum_exp(wvec);
  }

  // Penalties --------------------------------------------------------------
  obj = obj + pen_bias(bias, mean_bias, var_bias);
  obj = obj + pen_seq_error(seq, mean_seq, var_seq);
  obj = obj + pen_seq_error(od, mean_od, var_od);
  return obj;
}


//' Objective function optimized by \code{\link{uni_em}}.
//'
//' @inheritParams uni_em
//' @param pivec The current parameters.
//'
//' @author David Gerard
//'
//' @return The objective optimized by \code{\link{uni_em}} during
//'     that separate unimodal EM algorithm.
//'
// [[Rcpp::export]]
double uni_obj(arma::vec pivec,
               arma::vec weight_vec,
               arma::mat lmat,
               long double lambda) {
  arma::vec lpi = lmat.t() * pivec;
  double obj = 0.0;
  for (int k = 0; k < weight_vec.n_elem; k++) {
    if ((weight_vec(k) > TOL) && (lpi(k) > TOL)) {
      obj = obj + weight_vec(k) * std::log(lpi(k));
    } else if ((weight_vec(k) > TOL * 1000) && (lpi(k) < TOL)) {
      obj = R_NegInf;
      break;
    } else {
      // do nothing.
    }
  }

  // Add the penalty ----
  double pen = 0.0;
  if (lambda > TOL) {
    pen = lambda * arma::sum(arma::log(pivec));
  }

  return obj + pen;
}

//' EM algorithm to fit weighted ash objective.
//'
//' Solves the following optimization problem
//' \deqn{\max_{\pi} \sum_k w_k \log(\sum_j \pi_j \ell_jk).}
//' It does this using a weighted EM algorithm.
//'
//' @param weight_vec A vector of weights. Each element of \code{weight_vec} corresponds
//'     to a column of \code{lmat}.
//' @param lmat A matrix of inner weights. The columns are the "individuals" and the rows are the "classes."
//' @param pi_init The initial values of \code{pivec}. Each element of \code{pi_init}
//'     corresponds to a row of \code{lmat}.
//' @param itermax The maximum number of EM iterations to take.
//' @param obj_tol The objective stopping criterion.
//' @param lambda The penalty on the pi's. Should be greater than 0 and really really small.
//'
//'
//' @return A vector of numerics. The update of \code{pivec} in
//'     \code{\link{flexdog_full}}.
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
arma::vec uni_em(arma::vec weight_vec,
                 arma::mat lmat,
                 arma::vec pi_init,
                 long double lambda,
                 int itermax,
                 double obj_tol) {
  // check input ----------------------------------------------
  int ploidy = weight_vec.n_elem - 1;
  if (obj_tol < TOL) {
    Rcpp::stop("uni_em: obj_tol should be greater than 0.");
  }
  if (itermax < 0) {
    Rcpp::stop("uni_em: itermax should be greater than or equal to 0.");
  }
  if (weight_vec.n_elem != pi_init.n_elem) {
    Rcpp::stop("uni_em: weight_vec and pi_init should have the same number of elements.");
  }
  if (lmat.n_rows != weight_vec.n_elem) {
    Rcpp::stop("uni_em: lmat should have ploidy + 1 rows.");
  }
  if (lmat.n_cols != weight_vec.n_elem) {
    Rcpp::stop("uni_em: lmat should have ploidy + 1 columns.");
  }
  if (lambda < 0.0) {
    Rcpp::stop("uni_em: lambda cannot be negative.");
  }

  // Run EM ---------------------------------------------------
  int index       = 0;
  double err      = obj_tol + 1.0;
  arma::vec pivec = pi_init;
  double obj      = uni_obj(pivec, weight_vec, lmat, lambda);
  double old_obj  = obj;
  double lsum     = 0.0;
  arma::mat etamat(ploidy + 1, ploidy + 1);
  arma::vec nvec(ploidy + 1);

  while ((index < itermax) & (err > obj_tol)) {
    old_obj = obj;
    // get eta_jk -----------------
    for (int k = 0; k <= ploidy; k++) {
      lsum = 0.0;
      for (int j = 0; j <= ploidy; j++) {
        etamat(j, k) = pivec(j) * lmat(j, k);
        lsum = lsum + etamat(j, k);
      }
      for (int j = 0; j <= ploidy; j++) {
        etamat(j, k) = etamat(j, k) / lsum;
      }
    }
    // get n_j's --------------------------------
    nvec = etamat * weight_vec + lambda;
    // normalize to get pi_j's ------------------
    pivec = nvec / arma::sum(nvec);
    // calculate objective and update stopping criteria
    obj = uni_obj(pivec, weight_vec, lmat, lambda);
    if (obj < old_obj - TOL) {
      break; // this almost always is OK because other bias starting points will work

      Rcpp::Rcout << "Index: "
                  << index
                  << std::endl
                  << "obj: "
                  << obj
                  << std::endl
                  << "pivec: "
                  << std::endl
                  << pivec.t()
                  << std::endl
                  << "weight_vec: "
                  << std::endl
                  << weight_vec.t()
                  << std::endl
                  << "lmat: "
                  << std::endl
                  << lmat
                  << std::endl
                  << "etamat: "
                  << std::endl
                  << etamat
                  << "nvec: "
                  << std::endl
                  << nvec.t()
                  << std::endl
                  << std::endl;
      Rcpp::stop("uni_em: Objective is not increasing.\n");
    }
    err = std::fabs(obj - old_obj);
    index++;

  }
  return pivec;
}


//' Objective for mixture of known dist and uniform dist.
//'
//' @param alpha The mixing weight.
//' @param pvec The known distribution (e.g. from assuming an
//'     F1 population).
//' @param weight_vec A vector of weights.
//'
//' @return The objective when updating \code{pivec} when \code{model = "f1"}
//'     or \code{model = "s1"} in \code{\link{flexdog_full}}.
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
double f1_obj(double alpha,
              arma::vec pvec,
              arma::vec weight_vec) {

  // check input ----------------------------------
  int ploidy = pvec.n_elem - 1;
  if (weight_vec.n_elem != ploidy + 1) {
    Rcpp::stop("f1_obj: pvec and weight_vec should be the same length.");
  }
  if ((alpha < 0) | (alpha > 1)) {
    Rcpp::stop("f1_obj: alpha should be between 0 and 1.");
  }

  // get obj --------------------------------------
  double obj = 0.0;
  double newp;
  for (int k = 0; k <= ploidy; k++) {
    newp = (1.0 - alpha) * pvec(k) + alpha / ((double)ploidy + 1.0);
    obj = obj + weight_vec(k) * std::log(newp);
  }
  return obj;
}
