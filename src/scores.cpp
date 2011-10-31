#include "scores.h"
//#include <RcppArmadillo.h>
//#include <Rcpp.h>

using namespace arma ;
using namespace Rcpp ;


void sampleScores(arma::mat& Z, arma::mat& A, arma::mat& F, Rcpp::NumericVector& sigma2inv, 
									int n, int k, int p) {
    arma::vec s2i(sigma2inv.begin(), p, false);
    arma::mat Sigma_inv = arma::diagmat(s2i);
    
    arma::mat U = arma::randn(k, n);
    arma::mat fcov = arma::inv(arma::eye(k,k) + arma::trans(A)*Sigma_inv*A);
    arma::mat R = arma::chol(fcov);
    F = R*U + fcov*arma::trans(A)*Sigma_inv*Z;
}


void sampleScoresK(arma::mat& Z, arma::mat& A, arma::mat& F, arma::colvec& tauinv_score, int n, int k) {
    arma::mat U = arma::randn(k, n);
    arma::mat fcov = arma::inv(diagmat(tauinv_score) + arma::trans(A)*A);
    arma::mat R = arma::chol(fcov);
    F = R*U + fcov*arma::trans(A)*Z;
}

void sampleScoresVar(arma::colvec& tauinv_score, arma::mat& F, double& taua, double& taub, int n, int k) {
	
	//arma::rowvec A_maxnnz_col = arma::sum(A_restrict, 0);
	arma::colvec sumsq = sum(square(F),1);
	double nd = double(n);
    for (int i=0; i<k; i++) {
        double tauafc = taua + 0.5*nd;
        double taubfc = taub + 0.5*sumsq(i);
        tauinv_score(i) = Rf_rgamma(tauafc, 1.0)/taubfc;
    }
}


SEXP updateScoresC( SEXP Z_, SEXP A_, SEXP F_ ){
    BEGIN_RCPP
    //fix to allow obs var neq 1
    Rcpp::NumericMatrix Fr(F_);    
    Rcpp::NumericMatrix Ar(A_); 
    Rcpp::NumericMatrix Zr(Z_);
    
    int n = Zr.ncol();
    int k = Ar.ncol();
    int p = Ar.nrow();
    
    Rcpp::NumericVector sigma2inv(p, 1.0);
    
    arma::mat Z(Zr.begin(), p, n, false); 
    arma::mat A(Ar.begin(), p, k, false);
    arma::mat F(Fr.begin(), k, n, false);
    
    sampleScores(Z, A, F, sigma2inv, n, k, p);
    
    END_RCPP
}

