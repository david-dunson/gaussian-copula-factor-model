#include "main_mcmc.h"
#include "loadings_sparse.h"
#include "loadings_gdp.h"
#include "erl.h"
#include "loadings.h"
#include "scores.h"
#include "imh.h"

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>

using namespace arma ;
using namespace Rcpp ;


SEXP MCMCstep( SEXP Z_, SEXP A_, SEXP F_, SEXP tauinv_, SEXP rho_ , 
        SEXP Ra_, SEXP maxes_, SEXP argsorts_, SEXP A_restrict_, 
        SEXP priors_, SEXP nsim_, SEXP nburn_, SEXP thin_, 
        SEXP printstatus_, SEXP keepscores_, SEXP keeploadings_, 
        SEXP more_args_){
  
	BEGIN_RCPP
	
	Rcpp::List prop_cov, prop_mean, prop_prec;
	Rcpp::NumericVector nlevels;
	
	GetRNGstate();
	
	// Set the Armadillo seed from R's 
	int seed = (int)Rf_runif(0.0, 10000.0);
	std::srand(seed);
	
	// Pull together extra arguments
	Rcpp::List arg(more_args_);
	
	bool quiet;
	try {
		quiet = as<bool>(arg["quiet"]);
	} catch(...) {
		quiet = false;
	}
	
	int method = as<int>(arg["method"]);
  int px  = as<int>(arg["px"]);
  int imh = as<int>(arg["imh"]);
	
	double loadings_v = as<double>(arg["loadings_var"]);
	
	int t_model = as<int>(arg["factor_scales"]);
	int save_max      = as<int>(arg["save_max"]);
	double df         = as<double>(arg["df"]);
	bool positive     = as<bool>(arg["positive"]);
	bool keepscores   = as<bool>(keepscores_);
	bool keeploadings = as<bool>(keeploadings_);
	int  printstatus  = as<int>(printstatus_);

	if (!quiet) { 
    Rcout << "Using loadings prior: ";
	  switch (method) {
	    case 0:
        Rcout << "normal\n"; break;
      case 1:
        Rcout << "pointmass\n"; break;
      case 2:
        Rcout << "gdp\n"; break;
    }
    Rprintf("Using PX: %s\n", (px > 0) ? "Yes" : "No"); 
	  Rprintf("Using loadings prior variance %4.2f\n", loadings_v); 
    Rprintf("Factor Scale Hyperparameters (MV t model): %s\n", (t_model > 0) ? "Yes" : "No"); 
    Rprintf("IMH: %s\n", (imh > 0) ? "Yes" : "No"); 
	}	

  //Indices of missing data, to be imputed
  SEXP mInd_ = arg["mInd"];
  IntegerMatrix mInd(mInd_);
  int nmis = mInd.nrow();

  //Setup Gaussian margins, if any
  SEXP evi = arg["error_var_i"];
	Rcpp::NumericVector error_var_i(evi);
	double numobs = std::accumulate(error_var_i.begin(), error_var_i.end(), 0.0);
	bool sample_error_var = true;
	if (numobs<1.0) {
		sample_error_var = false;
	}
	
  //Setup IMH for mixed FA
	if (imh>0) {
		prop_cov  = arg["prop.cov"];
		prop_mean = arg["prop.mean"];
		prop_prec = arg["prop.prec"];
		nlevels   = arg["nlevels"];
	}
	
	// Prior hyperparams
	Rcpp::List priors(priors_);
	double taua  = as<double>(priors["taua"]);
	double taub  = as<double>(priors["taub"]);
	double rhoa  = as<double>(priors["rhoa"]);
	double rhob  = as<double>(priors["rhob"]);
	double s2a   = as<double>(priors["s2a"]);
	double s2b   = as<double>(priors["s2b"]);
	double alpha = as<double>(priors["alpha"]);
  double beta  = as<double>(priors["beta"]);
  double bp    = as<double>(priors["bp"]);
  double bq    = as<double>(priors["bq"]);
	
	// Create Rcpp objects
	Rcpp::NumericMatrix A_restrictr(A_restrict_);
	Rcpp::IntegerMatrix Ra(Ra_);
	Rcpp::IntegerMatrix argsorts(argsorts_);
	Rcpp::NumericMatrix maxes(maxes_);
	Rcpp::NumericMatrix Fr(F_);
	Rcpp::NumericMatrix Ar(A_);
	Rcpp::NumericMatrix Zr(Z_);
	Rcpp::NumericVector tauinv(tauinv_);
	Rcpp::NumericVector rho(rho_);
	
	int n = Zr.ncol();
	int k = Ar.ncol();
	int p = Ar.nrow();
		
	Rcpp::NumericVector tauinv_px(p, 1.0);
	Rcpp::NumericVector sigma2inv(p, 1.0);
	Rcpp::NumericMatrix pnzR(p,k);
	
	// Make  arma objects
	arma::mat A_restrict(A_restrictr.begin(), p, k, false);
	arma::rowvec A_maxnnz_col = arma::sum(A_restrict, 0);
	arma::colvec A_maxnnz_row = arma::sum(A_restrict, 1);
	
	arma::mat Z(Zr.begin(), p, n, false);
	arma::mat A(Ar.begin(), p, k, false); 
	arma::mat F(Fr.begin(), k, n, false);
	arma::colvec tauinv_score(k);
	tauinv_score.fill(1.0);
	
	arma::mat A_prior_var(1,k);
	A_prior_var.fill(1.0);
	if (t_model==0) {
		A_prior_var.fill(loadings_v);
	}
	
	arma::mat Psi_inv;
	if (method==3) {
	  Psi_inv = pow(randu(p,k),-1.0);
	}
	
	// MCMC simulation parameters/containers
	int nsim  = as<int>(nsim_);
	int nburn = as<int>(nburn_);
	int thin  = as<int>(thin_);
	int samples = nsim/thin;
	arma::vec imh_acc(p);
	imh_acc.fill(0.0);

	// Output objects
	arma::colvec Aravel(Ar.begin(), p*k, false);
	arma::running_stat_vec<double> Apstat;
	
	arma::colvec Fravel(Fr.begin(), k*n, false);
	arma::running_stat_vec<double> Fpstat;
	
	arma::colvec pnz(pnzR.begin(), k*p, false);
	arma::running_stat_vec<double> pnzstat;
	
	int ii = maxes.ncol();
	int jj = maxes.nrow();
	arma::colvec maxravel(maxes.begin(), ii*jj, false);
	
	int isamp = 0;
	
	int max_size = (save_max > 0) ? samples*ii*jj : 1;
	Rcpp::NumericVector post_max(max_size);
	vec_iterator ipost_max = post_max.begin();
	
	int score_size = keepscores ? samples*k*n : 1;
	Rcpp::NumericVector Fp(score_size);
	vec_iterator iFp = Fp.begin();
	
	int load_size = keeploadings ? samples*k*p : 1;
	Rcpp::NumericVector Ap(load_size);
	vec_iterator iAp = Ap.begin();
	
	Rcpp::NumericMatrix sigma2inv_samp(samples, p);
	
	Rcpp::NumericVector scAp(k*p);
	vec_iterator iscAp = scAp.begin();
	arma::colvec armascAp(scAp.begin(), p*k, false);
	
	// Iterators
	mat_iterator iF = F.begin();
	mat_iterator iA = A.begin();
	mat_iterator ipnzR = pnzR.begin();
	mat_iterator imax = maxes.begin();
	
	//temporary matrix for applying reduction fcn
	arma::mat Atmp(p,k);
	Atmp.fill(0.0);
	arma::mat Ftmp(p,k);
	Ftmp.fill(0.0);
	
	mat_iterator iFtmp = Ftmp.begin();
	mat_iterator iAtmp = Atmp.begin();
		
////////////////////////////////////////////////////////////////////////////////////////////////
//
// Start MCMC
//	
////////////////////////////////////////////////////////////////////////////////////////////////
	if (!quiet) { Rprintf("Beginning MCMC...\n\n"); }
	for (int i=0; i<nburn+nsim; i++) {
		R_CheckUserInterrupt();
		
		switch(method) {		
			case 0: //Normal loadings
				sampleLoadings(Z, A, F, sigma2inv, error_var_i, A_restrict, loadings_v, px, n, p, k, positive);
				break;
				
			case 1: //Sparse loadings
				sampleSparseLoadingsJ(Z, A, F, tauinv_px, sigma2inv, error_var_i, rho, 
															A_restrict, pnzR, n, p, k, px );
				sampleRho(rho, Ar, A_restrict, rhoa, rhob);
				break;
				
			case 2: //GDP loadings
				sampleLoadingsGDP(Z, A, F, sigma2inv, error_var_i, A_restrict, loadings_v, alpha, beta,
												  px,  n,  p,  k, positive );
				break;
			case 3: //Beta smb
			  break;
		}
		
		if (t_model==0) {
			sampleScores(Z, A, F, sigma2inv, n, k, p) ;
		} else {
			sampleScoresK(Z, A, F, tauinv_score, n, k);
			sampleScoresVar(tauinv_score, F, taua, taub, n, k);
		}
		
		arma::mat ldata_mean = A*F;
		
		if (imh>0) {
      do_imh(Zr, A, F, Ra, maxes, argsorts, ldata_mean, sigma2inv,
                  df, prop_cov, prop_mean, prop_prec, nlevels, 
                  error_var_i,imh_acc, s2a, s2b);
	 	} else {
			if (!sample_error_var) {
				sampleZ(Zr, Ra, maxes, argsorts, ldata_mean);
			} else {
				for (int pp=0; pp<p; pp++) {
					if (error_var_i(pp) > 0.0) {
						//Normal marginal distribution; sample idiosyncratic variance
						double ssq = accu(square(Z.row(pp) - ldata_mean.row(pp)));
    				sigma2inv(pp) =  Rf_rgamma(0.5*n + s2a, 1.0)/(0.5*ssq + s2b);
					} else {
						//Update via ERL
						sampleZ(Zr, Ra, maxes, argsorts, ldata_mean, pp);
					}
				}
			}
		}
	  
    //impute Gaussian margins
    if(nmis>0) {
      for(int ii=0; ii<nmis; ii++) {
        int mi = mInd(ii,0)-1;
        int mj = mInd(ii,1)-1;
        Z(mi, mj) = ldata_mean(mi, mj) + sqrt(1./sigma2inv(mi))*Rf_rnorm(0,1);
      }
    }
    
			
	  if (i>1 && i%printstatus==0 &&!quiet) {
			Rprintf("iteration %d\n",i);
			if (imh>0) {
        Rprintf("IMH acceptance rates:\n");
				for(int pp=0; pp<p; pp++) {
					Rprintf("%.3f ", imh_acc[pp]/i);
				}
			Rprintf("\n");
			}
	  }
	
	if (i>=nburn && i%thin==0){
			
			//Rprintf("saving\n\n");
			//cacluate scaling matrix
			arma::colvec b = arma::sum(arma::pow(A,2), 1);
			
			//save scaled loadings
			arma::colvec btmp;
			if (t_model>0) {
			  Atmp = A*diagmat(pow(tauinv_score, -0.5));
			  Ftmp = diagmat(sqrt(tauinv_score))*F;
			  btmp = arma::sum(arma::pow(Atmp,2), 1);
			}
			
			if (keeploadings) {
					for (int j = 0; j<k*p; j++) {
						if (t_model==0) {
								iAp[isamp*k*p + j] = iA[j]; ///sqrt(b[j%p]);
								iscAp[j] = iA[j]/sqrt(1/sigma2inv[j%p] + b[j%p]);
							} else {
								iAp[isamp*k*p + j] = iAtmp[j]; ///sqrt(b[j%p]);
								iscAp[j] = iAtmp[j]/sqrt(1/sigma2inv[j%p] + btmp[j%p]);
							}
					}
			} else {
					for (int j = 0; j<k*p; j++) {
						if (t_model==0) {
								iscAp[j] = iA[j]/sqrt(1/sigma2inv[j%p]+b[j%p]);
							} else {
								iscAp[j] = iAtmp[j]/sqrt(1/sigma2inv[j%p] +btmp[j%p]);
							}
					}
			}
			
			
			if (keepscores) {
					for (int j = 0; j<k*n; j++) {
						if (t_model==0) {
								iFp[isamp*k*n + j] = iF[j];
							} else {
								iFp[isamp*k*n + j] = iFtmp[j];
							}
					}
			}
			
			if (save_max > 0) {
				for (int j=0; j<ii*jj; j++) {
					ipost_max[isamp*ii*jj + j] = imax[j];
				}
			}
			
			Fpstat(Fravel);
			Apstat(armascAp);
			pnzstat(pnz);
			
			Rcpp::NumericMatrix::Row srow = sigma2inv_samp(isamp, _);	
			srow = sigma2inv;
			
			isamp++;
		}
	}
	
	if (!keepscores)   { Fp[0]       = NA_REAL; }
	if (!keeploadings) { Ap[0]       = NA_REAL; }
	if (save_max<1)    { post_max[0] = NA_REAL; }
	
	return List::create(
					_["maxp"]    = post_max,
					_["sigma2inv"] = sigma2inv_samp,
					_["pnz"]     = pnzstat.mean(),
					_["Ap"]      = Ap,
					_["Ap.mean"] = Apstat.mean(),
					_["Ap.var"]  = Apstat.var(),
					_["Fp"]      = Fp,
					_["Fp.mean"] = Fpstat.mean(),
					_["Fp.var"]  = Fpstat.var());

	PutRNGstate();
	
	END_RCPP
}
