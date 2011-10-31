#ifndef _bfa_scores_h
#define _bfa_scores_h

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "bfa_common.h"

RcppExport SEXP updateScoresC( SEXP Z_, SEXP A_, SEXP F_ );

void sampleScores(arma::mat& Z, arma::mat& A, arma::mat& F,  
									Rcpp::NumericVector& sigma2inv, int n, int k, int p);
void sampleScoresVar(arma::colvec& tauinv_score, arma::mat& F, double& taua, double& taub, int n, int k);
void sampleScoresK(arma::mat& Z, arma::mat& A, arma::mat& F, arma::colvec& tauinv_score, int n, int k);

#endif
