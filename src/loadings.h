#ifndef _bfa_loadings_h
#define _bfa_loadings_h

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "bfa_common.h"

void sampleLoadings(arma::mat& Z, arma::mat& A, arma::mat& F,  
          Rcpp::NumericVector& sigma2inv,
					Rcpp::NumericVector& error_var_i,
					arma::mat& A_restrict, double A_prior_var, 
					int px, int n, int p, int k, bool positive=false,
					arma::mat prior_prec = arma::mat(1,1));

#endif
