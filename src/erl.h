#ifndef _bfa_erl_h
#define _bfa_erl_h

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "bfa_common.h"

void sampleZ( Rcpp::NumericMatrix& Z , Rcpp::IntegerMatrix& Ra, 
			 Rcpp::NumericMatrix& maxes, Rcpp::IntegerMatrix& argsorts, 
			 arma::mat& mu, int row=-1);


RcppExport SEXP updateZ( SEXP Z_, SEXP Ra_, SEXP maxes_, SEXP argsorts_, SEXP A_, SEXP F_ );

#endif
