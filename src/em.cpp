#include "mupdog.h"

// [[Rcpp::export]]
void wem_fixp(arma::vec & pivec,
              const arma::mat & lmat,
              arma::mat & etamat,
              const int & nind,
              const int & nclass,
              const arma::vec & weight_vec,
              arma::vec & nvec,
              const long double & lambda) {
  double lsum = 0.0;
  // get eta_jk -----------------
  for (int k = 0; k < nind; k++) {
    lsum = 0.0;
    for (int j = 0; j < nclass; j++) {
      etamat(j, k) = pivec(j) * lmat(j, k);
      lsum = lsum + etamat(j, k);
    }
    for (int j = 0; j < nclass; j++) {
      etamat(j, k) = etamat(j, k) / lsum;
    }
  }
  // get n_j's --------------------------------
  nvec = etamat * weight_vec + lambda;
  // normalize to get pi_j's ------------------
  pivec = nvec / arma::sum(nvec);
}



// [[Rcpp::export]]
double wem_obj(const arma::vec & pivec,
               const arma::vec & weight_vec,
               const arma::mat & lmat,
               const long double & lambda) {
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


//' Generalized version of \code{\link{uni_em}}.
//'
//' @inherit uni_em
//'
//' @author David Gerard
//'
//' @export
//'
// [[Rcpp::export]]
arma::vec wem(arma::vec weight_vec,
              arma::mat lmat,
              arma::vec pi_init,
              long double lambda,
              int itermax,
              double obj_tol) {
  // check input ----------------------------------------------
  int nind   = weight_vec.n_elem;
  int nclass = pi_init.n_elem;
  if (obj_tol < TOL) {
    Rcpp::stop("uni_em: obj_tol should be greater than 0.");
  }
  if (itermax < 0) {
    Rcpp::stop("uni_em: itermax should be greater than or equal to 0.");
  }
  if (lmat.n_rows != pi_init.n_elem) {
    Rcpp::stop("uni_em: lmat should have ploidy + 1 rows.");
  }
  if (lmat.n_cols != weight_vec.n_elem) {
    Rcpp::stop("uni_em: lmat should have ploidy + 1 columns.");
  }
  if (lambda < 0.0) {
    Rcpp::stop("uni_em: lambda cannot be negative.");
  }

  // Run SQUAREM ------------------------------------------------
  int index       = 0;
  double err      = obj_tol + 1.0;
  arma::vec pivec = pi_init;
  double obj      = wem_obj(pivec, weight_vec, lmat, lambda);
  double old_obj  = obj;
  arma::mat etamat(nclass, nind);
  arma::vec nvec(nclass);
  arma::vec rn(nclass);
  arma::vec vn(nclass);
  double alpha_mpe1;
  double alpha_rre1;
  double alpha_hyb;
  double dot_rnvn;
  double wn;
  arma::vec pivec_n1(nclass);
  arma::vec pivec_n2(nclass);

  while ((index < itermax) & (err > obj_tol)) {
    old_obj = obj;
    pivec_n2 = pivec;

    // update pivec twice
    wem_fixp(pivec_n2,
             lmat,
             etamat,
             nind,
             nclass,
             weight_vec,
             nvec,
             lambda);

    pivec_n1 = pivec_n2; // one update

    wem_fixp(pivec_n2,
             lmat,
             etamat,
             nind,
             nclass,
             weight_vec,
             nvec,
             lambda);
    // both updates done.

    // get rn and vn ------------------------
    rn = pivec_n1 - pivec;
    vn = pivec_n2 - 2.0 * pivec_n1 + pivec;


    // get alphas ---------------------------
    dot_rnvn = arma::dot(rn, vn);
    alpha_mpe1 = arma::dot(rn, rn) / dot_rnvn;
    alpha_rre1 = dot_rnvn / arma::dot(vn, vn);
    wn = std::sqrt(alpha_rre1 / alpha_mpe1);
    alpha_hyb = wn * alpha_mpe1 + (1.0 - wn) * alpha_rre1;

    // get update ---------------------------
    pivec = pivec - (2.0 * alpha_hyb) * rn + std::pow(alpha_hyb, 2.0) * vn;

    for (int k = 0; k < nclass; k++) {
      if (pivec(k) < 0.0) {
        pivec(k) = 0.0;
      }
    }

    // one last fixed point -----------------
    wem_fixp(pivec,
             lmat,
             etamat,
             nind,
             nclass,
             weight_vec,
             nvec,
             lambda);

    // calculate objective and update stopping criteria
    obj = wem_obj(pivec, weight_vec, lmat, lambda);

    err = std::abs(obj - old_obj);
    index++;
  }
  return pivec;
}
