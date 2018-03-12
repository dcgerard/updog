#include "mupdog.h"

// functions for fitting flexdog

// [[Rcpp::export]]
NumericVector get_probk_vec(NumericVector pivec, std::string model, double mode) {
  int K = pivec.length() - 1;
  NumericVector probk_vec(K + 1);
  if (model == "flex") {
    probk_vec = pivec;
  } else if (model == "ash") {
    double denom; // what you divide the pi's by.
    for (int i = 0; i <= K; i++) { // iterate through pivec
      if (abs((double)i - mode) < TOL) {
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
    Rcpp::stop("get_probk_vec: `model` must be one of `'flex'` or `'ash'`");
  }
  return probk_vec;
}
