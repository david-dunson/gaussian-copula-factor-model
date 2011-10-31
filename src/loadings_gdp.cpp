#include "loadings_gdp.h"
#include "loadings.h"
#include "bfa_common.h"
//#include <RcppArmadillo.h>
//#include <Rcpp.h>

using namespace arma ;
using namespace Rcpp ;


void sampleLoadingsGDP(arma::mat& Z, arma::mat& A, 
          arma::mat& F, Rcpp::NumericVector& sigma2inv,
					Rcpp::NumericVector& error_var_i,
					arma::mat& A_restrict, double A_prior_var, 
					double alpha, double beta, int px, int n, int p, int k, bool positive ){
    
    arma::mat Lam = arma::zeros(p,k);
    //sample lambda, top level mixing param
    for (int i=0; i<p; i++) {
      for (int j=0; j<k; j++) {
        if (A_restrict(i,j)==0) {
          Lam(i,j) = 0.0;
        } else {
          Lam(i,j) = Rf_rgamma(alpha + 1, 1.0/(beta + fabs(A(i,j))));
        }
      }
    }
    
    arma::mat Psi_inv = arma::zeros(p,k);
    //sample Psi inverse, lower level mixing param
    for (int i=0; i<p; i++) {
      for (int j=0; j<k; j++) {
        if (A_restrict(i,j)==0) {
          Psi_inv(i,j) = 0.0;
        } else if (A_restrict(i,j)==2) {
          Psi_inv(i,j) = 1.0/A_prior_var; //truncated normal, not GDP.
        } else {
          double mu = fabs(A(i,j)*sqrt(Lam(i,j)));
          Psi_inv(i,j) = rinvgauss(mu, Lam(i,j)*Lam(i,j));
        }
      }
    }
    
    sampleLoadings(Z, A, F, sigma2inv, error_var_i, A_restrict, 1.0, px, 
                   n, p, k, positive, Psi_inv);
}

