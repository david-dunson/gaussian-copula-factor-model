#include "erl.h"
//#include <RcppArmadillo.h>
//#include <Rcpp.h>

using namespace arma ;
using namespace Rcpp ;

void sampleZ( Rcpp::NumericMatrix& Z , Rcpp::IntegerMatrix& Ra, Rcpp::NumericMatrix& maxes, 
							Rcpp::IntegerMatrix& argsorts, arma::mat& mu, int row){
	
	int cn;
	
	int n = Ra.ncol(), p = Ra.nrow();
	
	double etol = 1e-13;
	double cur_max, cur_min, next_min;
	int k, first_rank, cur_rank;
	double cm = -7.5;
	
	int slice = 0;
	int iter = 3;
	
	for (int i = 0; i < p; i++) {
		if (row==i || row < 0) {
			//double s = vars[i];
			first_rank = 1;
			cur_min  = -70.7;//*sqrt(s);
			next_min = -70.7;//*sqrt(s);
			cur_max = maxes(i, 0);
			cur_rank = 0;
			for (int j = 0; j < n; j++) {
				k = argsorts(i,j);
				if (cur_min>cur_max) {
					//Rprintf("bounds check failed @ row %d col %d: min = %f, max = %f, rank=%d\n\n", 
					//				i, k, cur_min, cur_max, Ra(i,k));
					Z(i,k) = cur_min;
				} else {
					if (Ra(i,k) < 0) {
					
						Z(i, k) = Rf_rnorm(mu(i,k), 1.0);
						cn = 1;
							
					}  else if (Ra(i,k)==0 && first_rank==1) {
					
						Z(i, k) = rtnorm_1(mu(i,k), 1.0, cur_min, cur_max);
						first_rank = 0;
						next_min = Z(i,k);
						cn = 2;
						
					} else if (cur_rank==Ra(i,k) and cur_rank==0){
					
						Z(i, k) = rtnorm_1(mu(i,k), 1.0, cur_min, cur_max);
						next_min = std::max(Z(i,k), next_min);
						cn = 3;
						
					} else if (cur_rank==Ra(i,k)) {
					
						Z(i, k) = rtnorm_1(mu(i,k), 1.0, cur_min, cur_max);
						next_min = std::max(Z(i,k), next_min);
						maxes(i, (cur_rank-1) ) = std::min(Z(i,k), maxes(i, (cur_rank-1) ));
						cn = 4;
					} else {
						cur_rank += 1;
						cur_min = next_min;
						cur_max = maxes(i, cur_rank);
						if (cur_max== 99999.0) {
							cur_max = 70.7;//R_PosInf;//*sqrt(s);
						}
						Z(i, k) = rtnorm_1(mu(i,k), 1.0, cur_min, cur_max);
						next_min = std::max(Z(i,k), next_min);
						maxes(i, (cur_rank-1) ) = Z(i,k);
						cn = 5;
					}
				}
			}
		}
	}
}


SEXP updateZ( SEXP Z_, SEXP Ra_, SEXP maxes_, SEXP argsorts_, SEXP A_, SEXP F_ ){
    BEGIN_RCPP
    GetRNGstate();
    
    Rcpp::IntegerMatrix Ra(Ra_);
    Rcpp::IntegerMatrix argsorts(argsorts_);
    Rcpp::NumericMatrix maxes(maxes_);
    Rcpp::NumericMatrix Fr(F_);
    Rcpp::NumericMatrix Ar(A_);
    Rcpp::NumericMatrix Zr(Z_);
    
    int n = Zr.ncol();
    int k = Ar.ncol();
    int p = Ar.nrow();
    
    arma::mat Z(Zr.begin(), p, n, false);      
    arma::mat A(Ar.begin(), p, k, false);      
    arma::mat F(Fr.begin(), k, n, false);       // reuses memory and avoids extra copy
    arma::mat mu = A*F;
    arma::vec vars = 1 + arma::sum(arma::pow(A,2), 1);
    sampleZ( Zr , Ra, maxes, argsorts, mu);
  
    PutRNGstate();
    END_RCPP
}

