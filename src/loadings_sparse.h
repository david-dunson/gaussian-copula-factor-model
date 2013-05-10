#ifndef _bfa_loadings_sparse_h
#define _bfa_loadings_sparse_h

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "bfa_common.h"

RcppExport SEXP updateSparseLoadingsJ( SEXP Z_, SEXP A_, SEXP F_, SEXP tauinv_, SEXP rho_, SEXP A_restrict_, SEXP pnz_ ) ;
RcppExport SEXP updateRho( SEXP rho_, SEXP A_, SEXP A_restrict_, SEXP rhoa_, SEXP rhob_ );

void sampleRho(Rcpp::NumericVector& rho, Rcpp::NumericMatrix& A, arma::mat A_restrict, 
			   double rhoa, double rhob);

void sampleSparseLoadingsJ(arma::mat& Z, arma::mat& A, arma::mat& F, 
						   Rcpp::NumericVector& tauinv, Rcpp::NumericVector& sigma2inv, Rcpp::NumericVector& error_var_i, 
						   Rcpp::NumericVector& rho, 
						   arma::mat& A_restrict, Rcpp::NumericMatrix& pnz, 
						   int n, int p, int k, int px, double A_prior_var=1.0 );

#endif
