#include "loadings_sparse.h"
#include "bfa_common.h"
//#include <RcppArmadillo.h>
//#include <Rcpp.h>

using namespace arma ;
using namespace Rcpp ;

void sampleRho(Rcpp::NumericVector& rho, Rcpp::NumericMatrix& A, arma::mat A_restrict, double rhoa, double rhob){
    int k = A.ncol();
    int p = A.nrow();
	arma::rowvec A_maxnnz_col = arma::sum(A_restrict, 0);
    for (int i=0; i<k; i++) {
        NumericMatrix::Column col = A.column(i) ;
        NumericVector As = ifelse(abs(col)<1e-10, 0.0, 1.0);
        double nnz = std::accumulate(As.begin(), As.end(), 0.0);
        double maxnnz = A_maxnnz_col(i);
        rho(i) = Rf_rbeta(rhoa + nnz, rhob + maxnnz-nnz);
    }
}

SEXP updateRho( SEXP rho_, SEXP A_, SEXP A_restrict_, SEXP rhoa_, SEXP rhob_ ){
    BEGIN_RCPP
    
    Rcpp::NumericVector rho(rho_);
    Rcpp::NumericMatrix A(A_);
    double rhoa = as<double>(rhoa_);
    double rhob = as<double>(rhob_);
	int k = A.ncol();
    int p = A.nrow();
	
	Rcpp::NumericMatrix A_restrictr(A_restrict_);
	arma::mat A_restrict(A_restrictr.begin(), p, k, false);
    
    sampleRho(rho, A, A_restrict, rhoa, rhob);
    
    END_RCPP
}


void sampleSparseLoadingsJ(arma::mat& Z, arma::mat& A, arma::mat& F, Rcpp::NumericVector& tauinv, 
							 Rcpp::NumericVector& sigma2inv, Rcpp::NumericVector& error_var_i,
						   Rcpp::NumericVector& rho, arma::mat& A_restrict, Rcpp::NumericMatrix& pnz, 
						   int n, int p, int k, int px, double A_prior_var ){
    
  arma::mat t, uvector, Ap_var, Ap_means, FFt;
  arma::vec sumf2 = arma::sum(arma::pow(F, 2), 1);
  double r = 1.0, samp = 0.0;

  if (px>0) {
	  FFt = F*arma::trans(F);
	  /*Ap_var = arma::inv(FFt + (1/A_prior_var)*arma::eye(k,k));
	  Ap_means = Ap_var*F*trans(Z);*/
  }
  
  
  double u=0.0, uu=0.0, v=0.0, loglog=0.0, logc=0.0;
  for (int i=0; i<p; i++) {
    
    //PX step
    if (px>0 && error_var_i(i)<1.0) {
      Ap_var = inv(FFt*sigma2inv(i) + (1/A_prior_var)*eye(k,k));
      Ap_means = Ap_var*F*trans(Z.row(i))*sigma2inv(i);
      double ssq = accu(square(trans(Z.row(i)) - trans(F)*Ap_means));
      double prior_ssq = as_scalar(trans(Ap_means) * Ap_means)/A_prior_var;
      double scale = 2.0/(ssq + prior_ssq);
      r = Rf_rgamma(0.5*n, scale);
    } else {
      r = 1.0;
    }
    
		for (int j=0; j<k; j++){
			if (A_restrict(i,j) > 0) { //nonzero
				t = Z.row(i) - A.row(i)*F + A(i,j)*F.row(j);
				u = arma::accu(F.row(j)%t)*sqrt(r)*sigma2inv(j);
				v = sumf2(j)*sigma2inv(j) + tauinv[i];
				if (A_restrict(i,j) > 1) { //strictly positive
					double lo = Rcpp::stats::pnorm_1( -u/sqrt(v), 0.0, true, false);
					double un = std::min(Rf_runif(lo, 1.0), 1.0-6.7e-16);
					double samp = u/v + Rcpp::stats::qnorm_1(un, 0.0, true, false)/sqrt(v);
					A(i,j) = samp;
				} else {
					loglog = log(u*u) - log(2.0*v);
					logc = log(rho[j]) - log(1-rho[j]) - 0.5*log(v) + 0.5*log(tauinv[i]) + exp(loglog);
					pnz(i,j) = 1/(1+exp(logc));
					uu = Rf_runif(0.0, 1.0);
					log(1/uu-1) > logc ? A(i,j) = 0.0 : A(i,j) = Rf_rnorm(u/v, 1/sqrt(v));
				}
			} else {
				A(i,j) = 0.0;
				pnz(i,j) = 1.0;
			}
		}
	}
}

SEXP updateSparseLoadingsJ( SEXP Z_, SEXP A_, SEXP F_, SEXP tauinv_, SEXP rho_, SEXP A_restrict_, SEXP pnz_ ){
    BEGIN_RCPP
    
    Rcpp::NumericMatrix Fr(F_);          
    Rcpp::NumericMatrix Ar(A_);           
    Rcpp::NumericMatrix Zr(Z_);            
    Rcpp::NumericVector tauinv(tauinv_);
    Rcpp::NumericVector rho(rho_);
    Rcpp::NumericMatrix pnz(pnz_);
	
		int n = Zr.ncol();
    int k = Ar.ncol();
    int p = Ar.nrow();
	
		Rcpp::NumericMatrix A_restrictr(A_restrict_);
		arma::mat A_restrict(A_restrictr.begin(), p, k, false);
		arma::rowvec A_maxnnz_col = arma::sum(A_restrict, 0);
		arma::colvec A_maxnnz_row = arma::sum(A_restrict, 1);
			
    GetRNGstate();
    
    arma::mat Z(Zr.begin(), p, n, false);       // reuses memory and avoids extra copy
    arma::mat A(Ar.begin(), p, k, false);       // reuses memory and avoids extra copy
    arma::mat F(Fr.begin(), k, n, false);       // reuses memory and avoids extra copy
    //fix to allow px
    //fix to allow obs var neq 1
    Rcpp::NumericVector sigma2inv(p, 1.0);
    Rcpp::NumericVector error_var_i(p, 0.0); 
    sampleSparseLoadingsJ(Z, A, F, tauinv, sigma2inv, error_var_i, rho, A_restrict, pnz, n, p, k, 0 );
    
    PutRNGstate();
    
    END_RCPP
}


