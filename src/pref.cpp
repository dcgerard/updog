#include "mupdog.h"
#include <sstream>

// functions to handle preferential pairing model

// gets probability vector from pairing vector.
// [[Rcpp::export]]
NumericVector dist_from_p(IntegerVector p, int ploidy) {
  if (p.length() != 3) {
    Rcpp::stop("dist_from_p: p.length() should be 3.");
  }
  if (Rcpp::sum(p) != (ploidy / 2)) {
    Rcpp::stop("dist_from_p: sum(p) should be ploidy / 2.");
  }

  NumericVector probvec(ploidy / 2 + 1);
  for (int ell = 0; ell <= p(1); ell++) {
    probvec(ell + p(2)) = R::dbinom((double)ell, (double)p(1), 0.5, false);
  }

  return probvec;
}

// [[Rcpp::export]]
int num_pairs(int ell, int ploidy) {
  int np;
  if (2 * ell >= ploidy) {
    np = ploidy / 2 - std::ceil((double)ell / 2) + 1;
  } else {
    np = std::floor((double)ell / 2) + 1;
  }
  return np;
}

// [[Rcpp::export]]
std::string get_pname(IntegerVector pvec) {
  std::ostringstream strs;
  strs << "("
       << pvec(0)
       << ","
       << pvec(1)
       << ","
       << pvec(2)
       << ")";
  std::string pname = strs.str();
  return pname;
}


//' Returns segregation probabilities, pairing representation
//' and number of ref alleles given the ploidy.
//'
//'
//' @param ploidy The ploidy of the individual. Should be even
//'     and greater than 0.
//'
//' @return A list of three elements
//' \describe{
//' \item{\code{probmat}}{The rows index the pairing configuration
//'     and the columns index the number of reference alleles
//'     segregating. The elements are the probability of
//'     segregating the given number of reference alleles
//'     in a given category.}
//' \item{\code{pmat}}{The pairing representation of the
//'     configuration.}
//' \item{\code{lvec}}{The number of reference alleles an individual
//'     has given their pairing configuration in \code{pmat}.}
//' }
//'
//' @author David Gerard
//'
//' @export
//'
//' @examples
//' get_bivalent_probs(4)
//'
// [[Rcpp::export]]
List get_bivalent_probs(int ploidy) {
  if (ploidy <= 0) {
    Rcpp::stop("get_bivalent_probs: ploidy should be greater than 0.");
  }
  if ((ploidy % 2) == 1) {
    Rcpp::stop("get_bivalent_probs: ploidy needs to be even.");
  }

  int nrow = (ploidy + 4) * (ploidy + 2) / 8;
  int ncol = ploidy / 2 + 1;

  NumericMatrix probmat(nrow, ncol);
  IntegerMatrix pmat(nrow, 3);
  IntegerVector lvec(nrow);
  CharacterVector rowname_vec(nrow);

  int i = 0;
  IntegerVector pvec(3);
  int numcat = 0;
  for (int ell = 0; ell <= ploidy; ell++) {
    pvec(1) = ell % 2;
    pvec(2) = std::floor((double)ell / 2);
    pvec(0) = ploidy / 2 - pvec(1) - pvec(2);
    numcat = num_pairs(ell, ploidy);
    for (int j = 0; j < numcat; j++) {
      probmat(i, _) = dist_from_p(pvec, ploidy);
      rowname_vec(i) = get_pname(pvec);
      pmat(i, 0) = pvec(0);
      pmat(i, 1) = pvec(1);
      pmat(i, 2) = pvec(2);
      lvec(i) = ell;
      pvec(0) = pvec(0) - 1;
      pvec(1) = pvec(1) + 2;
      pvec(2) = pvec(2) - 1;
      i++;
    }
  }

  CharacterVector colname_vec(ncol);
  for (int i = 0; i < ncol; i++) {
    colname_vec(i) = (char)i;
  }

  rownames(probmat) = rowname_vec;
  colnames(probmat) = colname_vec;

  CharacterVector pcolname(3);
  pcolname(0) = "aa";
  pcolname(1) = "Aa";
  pcolname(2) = "AA";
  colnames(pmat) = pcolname;

  return List::create(Named("probmat") = probmat,
                      Named("pmat") = pmat,
                      Named("lvec") = lvec);
}



// counts the number of pairings where the chromosome labels matter
// but the pairing label does not matter.
// [[Rcpp::export]]
int count_pairings(int ploidy) {
  if (ploidy <= 0) {
    Rcpp::stop("count_pairings: ploidy should be greater than 0.");
  }
  if ((ploidy % 2) == 1) {
    Rcpp::stop("count_pairings: ploidy should be even.");
  }

  double lratio = 0.0;
  for (int i = 2; i <= ploidy; i += 2) {
    lratio += R::lchoose((double)i, 2.0);
  }
  lratio = lratio - R::lgammafn((double)ploidy / 2.0 + 1.0);

  return std::round(std::exp(lratio));
}

// Counts the number of ways to choose AA into p2 bins given ell A's.
// You can get the counts for aa into p0 bins given ploidy - ell a's
// by setting ell = ploidy - ell and p2 = p0.
// [[Rcpp::export]]
int count_doubles(int ell, int p2) {
  if (ell < 0) {
    Rcpp::stop("count_doubles: ell should be >= 0.");
  }
  if (p2 < 0) {
    Rcpp::stop("count_doubles: p2 should be >= 0.");
  }

  double lratio = 0.0;
  if (p2 == 0) {
    lratio = 0.0;
  } else if (2 * p2 > ell) {
    lratio = R_NegInf;
  } else {
    for (int i = ell - 2 * p2 + 2; i <= ell; i += 2) {
      lratio += R::lchoose((double)i, 2.0);
    }
    lratio = lratio - R::lgammafn((double)p2 + 1.0);
  }

  return std::round(std::exp(lratio));
}

// [[Rcpp::export]]
int count_pairings_given_p(int ploidy, IntegerVector pvec) {

  if (pvec.length() != 3) {
    Rcpp::stop("count_pairings_given_p: pvec should have length 3");
  }
  if (ploidy != 2 * Rcpp::sum(pvec)) {
    Rcpp::stop("count_pairings_given_p: ploidy should equal 2 * sum(pvec).");
  }
  if (ploidy <= 0) {
    Rcpp::stop("count_pairings_given_p: ploidy should be > 0.");
  }
  if ((ploidy % 2) == 1) {
    Rcpp::stop("count_pairings_given_p: ploidy should be even.");
  }

  int ell = pvec(1) + 2 * pvec(2);

  int x1 = count_doubles(ell, pvec(2));
  int x2 = count_doubles(ploidy - ell, pvec(0));
  int x3 = R::gammafn((double)pvec(1) + 1.0);

  int numpair = x1 * x2 * x3;
  return numpair;
}

//' Return mixture weights needed to obtain a hypergeometric
//' distribution.
//'
//' Obtains the mixing weights for the mixing distributions
//' of \code{\link{get_bivalent_probs}} to return a hypergeometric
//' distribution where \code{ploidy} is the population size,
//' \code{ell} is the number of success states in the population,
//' and \code{ploidy / 2} is the number of draws. If these
//' are the mixing weights in the population, then there is no
//' preferential pairing.
//'
//' @param ploidy The ploidy of the individual.
//' @param ell The number of reference alleles in the individual.
//'
//' @return A list with the following two elements:
//' \describe{
//' \item{\code{pmat}}{Reach row is a category and the columns
//'     index either aa, Aa, or AA.}
//' \item{\code{weightvec}}{The mixing weights for each row of pmat.}
//' }
//'
//' @author David Gerard
//'
//' @export
//'
//' @examples
//' get_hyper_weights(4, 2)
//'
// [[Rcpp::export]]
List get_hyper_weights(int ploidy, int ell) {
  if (ploidy <= 0) {
    Rcpp::stop("get_hyper_weights: ploidy should be > 0");
  }
  if ((ploidy % 2) == 1) {
    Rcpp::stop("get_hyper_weights: ploidy should be even");
  }
  if ((ell < 0) | (ell > ploidy)) {
    Rcpp::stop("get_hyper_weights: ell should be between 0 and ploidy (inclusive).");
  }

  if (ploidy > 20) {
    Rcpp::stop("get_hyper_weights: we run into numerical stability issues that we have not fixed yet when ploidy > 20.");
  }

  int np = count_pairings(ploidy); // number of pairings
  int np_p; //number of pairings given pvec.
  IntegerVector pvec(3);
  pvec(1) = ell % 2;
  pvec(2) = std::floor((double)ell / 2);
  pvec(0) = ploidy / 2 - pvec(1) - pvec(2);
  int numcat = num_pairs(ell, ploidy);
  IntegerMatrix pmat(numcat, 3);
  NumericVector mixing_weights(numcat);
  for (int j = 0; j < numcat; j++) {
    np_p = count_pairings_given_p(ploidy, pvec);
    mixing_weights(j) = (double)np_p / (double)np;
    pmat(j, 0) = pvec(0);
    pmat(j, 1) = pvec(1);
    pmat(j, 2) = pvec(2);
    pvec(0) = pvec(0) - 1;
    pvec(1) = pvec(1) + 2;
    pvec(2) = pvec(2) - 1;
  }

  CharacterVector pcolname(3);
  pcolname(0) = "aa";
  pcolname(1) = "Aa";
  pcolname(2) = "AA";
  colnames(pmat) = pcolname;

  return List::create(Named("pmat") = pmat,
                      Named("weightvec") = mixing_weights);
}

