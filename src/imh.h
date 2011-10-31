#ifndef _bfa_imh_h
#define _bfa_imh_h

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "bfa_common.h"

void sampleZ_cut( Rcpp::NumericMatrix& Z , Rcpp::IntegerMatrix& Ra, 
			 Rcpp::NumericMatrix& maxes, Rcpp::IntegerMatrix& argsorts, 
			 arma::mat& mu, int row=-1);

double likrat( Rcpp::NumericMatrix& Z , Rcpp::IntegerMatrix& Ra, 
			   arma::vec& current, arma::vec& propose,
			   Rcpp::IntegerMatrix& argsorts, arma::mat& A, arma::mat& F, int row, int le);
			   
void do_imh(Rcpp::NumericMatrix& Zr, 
	          arma::mat& A,   
	          arma::mat& F,
            Rcpp::IntegerMatrix& Ra, 
            Rcpp::NumericMatrix& maxes, 
            Rcpp::IntegerMatrix& argsorts,  
            arma::mat& ldata_mean,
            Rcpp::NumericVector& sigma2inv,
            double df,
            Rcpp::List& prop_cov, 
            Rcpp::List& prop_mean, 
            Rcpp::List& prop_prec,
            Rcpp::NumericVector& nlevels,
            Rcpp::NumericVector& error_var_i,
            arma::vec& imh_acc,
            double s2a,
            double s2b);
            
//SEXP updateZ( SEXP Z_, SEXP Ra_, SEXP maxes_, SEXP argsorts_, SEXP A_, SEXP F_ );
RcppExport SEXP updateZcut( SEXP Z_, SEXP Ra_, SEXP maxes_, SEXP argsorts_, SEXP A_, SEXP F_ );

#endif
