#include "imh.h"
//#include <RcppArmadillo.h>
//#include <Rcpp.h>

using namespace arma ;
using namespace Rcpp ;

void sampleZ_cut( Rcpp::NumericMatrix& Z , Rcpp::IntegerMatrix& Ra, 
                  Rcpp::NumericMatrix& maxes, Rcpp::IntegerMatrix& argsorts, 
                  arma::mat& mu, int row){

	int n = Ra.ncol(), p = Ra.nrow();
	
	double etol = 1e-13;
	double cur_max, cur_min, next_min;
	int k, first_rank, cur_rank;
	double cm = -7.5;
	
	int slice = 0;
	int iter = 3;
	
	for (int i = 0; i < p; i++) {
		if (row==i || row<0) {
			//double s = vars[i];
			first_rank = 1;
			cur_min  = -70.0;//R_NegInf; //-7.7*sqrt(s);
			next_min = -70.0;//R_NegInf; //-7.7*sqrt(s);
			cur_max = maxes(i, 0);
			cur_rank = 0;
			for (int j = 0; j < n; j++) {
				k = argsorts(i,j);
				if (cur_min>cur_max) {
					Rprintf("bounds check failed @ row %d col %d: min = %f, max = %f\n\n", i, k, cur_min, cur_max);
				}
				if (Ra(i,k) < 0) {
						Z(i, k) = Rf_rnorm(mu(i,k), 1.0);
				} else if (Ra(i,k)==0 && first_rank==1) {
					Z(i, k) = rtnorm_1(mu(i,k), 1.0, cur_min, cur_max);
					first_rank = 0;
				} else if (cur_rank==Ra(i,k) and cur_rank==0){
					Z(i, k) = rtnorm_1(mu(i,k), 1.0, cur_min, cur_max);
				} else if (cur_rank==Ra(i,k)) {
					Z(i, k) = rtnorm_1(mu(i,k), 1.0, cur_min, cur_max);
				} else {
					cur_rank += 1;
					cur_min = cur_max;
					cur_max = maxes(i, cur_rank);
					if (cur_max== 99999.0) {
						cur_max = 70.0;//R_PosInf; //cur_max*sqrt(s);
					}
					Z(i, k) = rtnorm_1(mu(i,k), 1.0, cur_min, cur_max);
				}
			}
    }
  }
}


SEXP updateZcut( SEXP Z_, SEXP Ra_, SEXP maxes_, SEXP argsorts_, SEXP A_, SEXP F_ ){
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
    
    arma::mat Z(Zr.begin(), p, n, false);       // reuses memory and avoids extra copy
    arma::mat A(Ar.begin(), p, k, false);       // reuses memory and avoids extra copy
    arma::mat F(Fr.begin(), k, n, false);       // reuses memory and avoids extra copy
    arma::mat mu = A*F;
    
    sampleZ_cut( Zr , Ra, maxes, argsorts, mu);
  
    PutRNGstate();
    END_RCPP
}


double likrat( Rcpp::NumericMatrix& Z , Rcpp::IntegerMatrix& Ra, 
         arma::vec& current, arma::vec& propose,
         Rcpp::IntegerMatrix& argsorts, arma::mat& A, arma::mat& F, int row, int le){
    
    int cn;
    //arma::mat mu = A*F;
    int n = Ra.ncol(), p = Ra.nrow();
    
    int maxr = le;
    
    double etol = 1e-13;
    double cur_max, cur_min, next_min;
    double cm = -7.5;
    
    double oldlik = 0.0;
    double newlik = 0.0;
    double mu = 0.0;
    
    double ma = 0.0;
    double mi = 0.0;
    double propma = 0.0;
    double propmi = 0.0;
    int k;
    //Rprintf("%d\n",k);
		for (int j=0; j<n; j++) {
			k = Ra(row, j);
			//Rprintf("%d/%d\n", k, maxr);
			mu = arma::dot(A.row(row), F.col(j));
			if (k==0) {
				//lower = -inf
				ma = current(k);
				propma = propose(k);
				oldlik += Rcpp::stats::pnorm_1(ma, mu, true, true);
				newlik += Rcpp::stats::pnorm_1(propma, mu, true, true);
			} else if (k==maxr) {
				//upper = inf
				mi = current(k-1);
				propmi = propose(k-1);
				oldlik += log(1.0-Rcpp::stats::pnorm_1(mi, mu, true, false));
				newlik += log(1.0-Rcpp::stats::pnorm_1(propmi, mu, true, false));
			} else if (k>0) {
				//(a,b)
				ma = current(k);
				propma = propose(k);
				mi = current(k-1);
				propmi = propose(k-1);
				oldlik += log(Rcpp::stats::pnorm_1(ma, mu, true, false)
								- Rcpp::stats::pnorm_1(mi, mu, true, false));
				newlik += log(Rcpp::stats::pnorm_1(propma, mu, true, false)
								- Rcpp::stats::pnorm_1(propmi, mu, true, false)); 
			}
		}
  double out = newlik - oldlik;
  return(out);
}



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
            double s2b) {
  int p = A.n_rows;
  int n = F.n_cols;
  arma::mat Z(Zr.begin(), p, n, false);
  for (int pp=0; pp<p; pp++) {
    if (error_var_i(pp)<1) {
      int le = nlevels[pp];
      if (le > 1) {
        Rcpp::NumericVector mu_      = prop_mean[pp];
        Rcpp::NumericMatrix covmat_  = prop_cov[pp];
        Rcpp::NumericMatrix precmat_ = prop_prec[pp];
        Rcpp::NumericVector current_ = maxes(pp, _);
        
        arma::vec mu(mu_.begin(), le, false);
        arma::mat covmat(covmat_.begin(), le, le, false);
        arma::mat precmat(precmat_.begin(), le, le, false);
        
        arma::vec propose = mu + arma::trans(covmat)*my_randn(le)/sqrt(Rf_rchisq(df));
        int maxprop = 10000;
        
        bool good = false;
        for (int ii=0; ii<maxprop; ii++) {
          if (is_sorted(propose)) {
            good = true;
            break;
          } else {
            propose = mu + arma::trans(covmat)*my_randn(le)/sqrt(Rf_rchisq(df));
          }
        }
        
        if (!good) {
          Rprintf("proposal sampling failed\n");
        }
        
        arma::vec current(current_.begin(), le, false);
        
        double a = df + arma::as_scalar( arma::trans(current-mu)*precmat*(current-mu) );
        double b = df + arma::as_scalar( arma::trans(propose-mu)*precmat*(propose-mu) );
        double propa = 0.5*(le+df)*(log(a) - log(b));
  
        double lika = likrat(Zr, Ra, current, propose, argsorts, A, F, pp, le);
  
        if ( (log(Rf_runif(0.0, 1.0)) < propa+lika) && good) {
          imh_acc[pp] += 1.0;
          for (int cc=0; cc<le; cc++) {
            maxes(pp,cc) = propose(cc);
          }
          sampleZ_cut(Zr, Ra, maxes, argsorts, ldata_mean, pp);
        }
      } else {
        double mu = as<double>(prop_mean[pp]);
        double v = as<double>(prop_cov[pp]);
        double current = maxes(pp,0);
        double propose = Rf_rnorm(mu, sqrt(v/Rf_rchisq(df)));
        double a = df + (current-mu)*(current-mu)/sqrt(v);
        double b = df + (propose-mu)*(propose-mu)/sqrt(v);
        double propa = 0.5*(le+df)*(log(a) - log(b));
        
        arma::vec acurrent(1);
        acurrent.fill(current);
        arma::vec apropose(1);
        apropose.fill(propose);
        double lika = likrat(Zr, Ra, acurrent, apropose, argsorts, A, F, pp, le);
        
        if (log(Rf_runif(0.0, 1.0)) < propa+lika) {
          imh_acc[pp] += 1.0;
          maxes(pp,0) = propose;
          sampleZ_cut(Zr, Ra, maxes, argsorts, ldata_mean, pp);
        }
      }
    } else {
      //Normal marginal distribution; sample idiosyncratic variance
      double ssq = accu(square(Z.row(pp) - ldata_mean.row(pp)));
      sigma2inv(pp) =  Rf_rgamma(0.5*n + s2a, 1.0)/(0.5*ssq + s2b);
    }
  }
}






