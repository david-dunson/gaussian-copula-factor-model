#include "loadings_beta.h"
#include "loadings.h"
#include "bfa_common.h"
//#include <RcppArmadillo.h>
//#include <Rcpp.h>

using namespace arma ;
using namespace Rcpp ;

double g(double psi, double a, double p, double q) {
  return -0.5*a*a/(psi*psi)+(1/p-2)*log(psi)+(q-1)*log(1-psi);
}

double sl(double ns, double psi0, double a, double p, double q) {
  int ii = 100;
  double xm = fabs(a)/sqrt(2-1/p);
  double yrt = -a*a*0.5;
  double step = 0.01; //std::min(0.01, 0.5*xm+0.001); //2nd term is sd from laplace approx
  
  double x0 = psi0;
  double y0 = g(x0, a, p, q);
  double x = x0;
  double y = y0;
  double gx, lo, hi, glo, ghi, xs, gxs;
  for (int i=0; i<ns; i++) {
    gx = g(x, a, p, q);
    y = gx - Rf_rexp(1);
    lo = x;
    hi = x;
    glo = y + 1.0;
    ghi = y + 1.0;
    // Get looooooooowwwwwww
    for (double s=0; s<ii; i++) {
      if (glo>y && fabs(lo)>1e-8) {
        lo = std::max(lo - (s+1)*step, 1e-8);
        glo = g(lo, a, p, q);
      } else {
        break;
      }
    }
    
    if (gx<yrt) {
      hi = 1.0;
    } else {
      for (double s=0; s<ii; s++) {
        if (ghi>y) {
          hi = std::min(hi+(s+1)*step, 1.0-1e-8);
          ghi = g(hi, a, p, q);
        } else {
          break;
        }
      }
    }
    xs = Rf_runif(lo,hi);
    gxs = g(xs, a, p, q);
    int maxi = 100000;
    int it = 0;
    while (gxs < y && it<maxi) {
      if (xs>x) { hi = xs; }
      if (xs<x) { lo = xs; }
      xs = Rf_runif(lo, hi);
      gxs = g(xs, a, p, q);
      it++;
      R_CheckUserInterrupt();
    }
    x = xs;
  }
  return x;
}

double betaslice(double psi0, double a, double p, double q) {
  double ns;
  if (fabs(a)<0.01 || fabs(a)>5) {
    ns=20.0;
  } else {
    ns=5.0;
  }
  return sl(ns, psi0, a, p, q);
}

void sampleLoadingsBeta(arma::mat& Z, arma::mat& A, 
          arma::mat& F, Rcpp::NumericVector& sigma2inv,
					Rcpp::NumericVector& error_var_i,
					arma::mat& A_restrict, double A_prior_var, 
					arma::mat &Psi_inv, double bp, double bq, 
					int px, int n, int p, int k, bool positive){
    
    double psi;
    //sample Psi inverse, lower level mixing param
    for (int i=0; i<p; i++) {
      for (int j=0; j<k; j++) {
        if (i < j) {
          Psi_inv(i,j) = 0.0;
        } else if (i==j) {
          Psi_inv(i,i) = 1.0/A_prior_var;
        } else {
          psi = betaslice(1.0/Psi_inv(i,j), A(i,j), bp, bq);
          Psi_inv(i,j) = 1.0/psi;
        }
      }
    }
    sampleLoadings(Z, A, F, sigma2inv, error_var_i, A_restrict, 1.0, px, 
                   n, p, k, positive, Psi_inv);
}

