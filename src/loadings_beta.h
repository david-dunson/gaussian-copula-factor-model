#ifndef _bfa_loadings_beta_h
#define _bfa_loadings_beta_h

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "bfa_common.h"

void sampleLoadingsBeta(arma::mat& Z, arma::mat& A, 
          arma::mat& F, Rcpp::NumericVector& sigma2inv,
					Rcpp::NumericVector& error_var_i,
					arma::mat& A_restrict, double A_prior_var, 
					arma::mat &Psi_inv, double bp, double bq, 
					int px, int n, int p, int k, 
					bool positive=false);

#endif
