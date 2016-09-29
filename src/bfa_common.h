#ifndef _bfa_common_h
#define _bfa_common_h

#include <RcppArmadillo.h>
#include <Rcpp.h>

typedef Rcpp::NumericMatrix::iterator mat_iterator;
typedef Rcpp::NumericVector::iterator vec_iterator;

bool is_sorted(arma::vec row);

bool bfa_isfinite(double x) ;
double rtnorm_slice(int iter, double mu, double sigma, double a, double b);
double rtnorm_1(double mu, double sigma, double a, double b);
double rinvgauss(double mu, double lam);

template<typename T, typename U> arma::mat my_randn(T k, U n) {
  Rcpp::NumericVector z = Rcpp::rnorm(k*n);
  arma::mat out(z.begin(), k, n);
  return out;
}
template<typename T, typename U> arma::mat my_randu(T k, U n) {
  Rcpp::NumericVector z = Rcpp::runif(k*n);
  arma::mat out(z.begin(), k, n);
  return out;
}

template<typename T> arma::mat my_randn(T k) {
  return my_randn(k, 1);
}
    
#endif
