#ifndef _bfa_common_h
#define _bfa_common_h

#include <RcppArmadillo.h>
#include <Rcpp.h>

typedef Rcpp::NumericMatrix::iterator mat_iterator;
typedef Rcpp::NumericVector::iterator vec_iterator;

bool is_sorted(arma::vec row);

bool isfinite(double x) ;
double rtnorm_slice(int iter, double mu, double sigma, double a, double b);
double rtnorm_1(double mu, double sigma, double a, double b);
double rinvgauss(double mu, double lam);

#endif
