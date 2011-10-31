#ifndef _bfa_mcmc_h
#define _bfa_mcmc_h

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "bfa_common.h"
//#include <google/profiler.h>

RcppExport SEXP MCMCstep( SEXP Z_, SEXP A_, SEXP F_, SEXP tauinv_, SEXP rho_ , 
			  SEXP Ra_, SEXP maxes_, SEXP argsorts_, SEXP A_restrict_, 
			  SEXP priors_, SEXP nsim_, SEXP nburn_, SEXP thin_, 
			  SEXP printstatus_, SEXP keepscores_, SEXP keeploadings_, 
			  SEXP more_args_);

#endif
