#include <RcppArmadillo.h>
#include <Rcpp.h>

#include "bfa_common.h"

using namespace arma ;
using namespace Rcpp ;

bool is_sorted(arma::vec row) {
	int len = row.n_elem;
	bool out = true;
	for (int i=1; i<len; i++) {
		if (row(i)-row(i-1)<0.0) {
			out=false;
			break;
		}
	}
	return(out);
}



bool isfinite(double x) 
{
	return (x <= DBL_MAX && x >= -DBL_MAX); 
}    



double rtnorm_slice(int iter, double mu, double sigma, double a, double b) {
	
	double up  = 0.0;
	double upp = 0.0;
	double lo = 0.0;
	double y  = 0.0;
	
	double x = b==R_PosInf ? a+Rf_runif(0, 1.0) : (b-a)/2; //(b-a)*Rf_runif(0.0, 1.0) + a;
	
	for (int i=0; i<iter; i++) {
		y = exp(-0.5*(x-mu)*(x-mu)/(sigma*sigma))*Rf_runif(0.0, 1.0);
		upp = mu + sqrt(-2.0*sigma*sigma*log(y));
		up =  (b==R_PosInf || upp<b) ? upp : b;
		lo = std::max(mu - sqrt(-2.0*sigma*sigma*log(y)), a);
		x = (up-lo)*Rf_runif(0.0, 1.0) + lo;
	}
	
	return(x);
	
}


double rtnorm_1(double mu, double sigma, double a, double b){
    //draw from mean mu, unit variance truncated normal
    //sigma argument is a farce
    double lo, hi, out;
    double zlo = a-mu; double zhi = b-mu;
    if(a<b) {
			if ((std::fabs(zlo) < 8.0) && (std::fabs(zhi) < 8.0)) {
				lo = Rcpp::stats::pnorm_1(a, mu, true, false);
				hi = Rcpp::stats::pnorm_1(b, mu, true, false);
			} else if (std::fabs(zlo) >8.0) {
				lo = 0.0;
				hi = Rcpp::stats::pnorm_1(b, mu, true, false);
			} else if (std::fabs(zhi) >8.0) {
				lo = Rcpp::stats::pnorm_1(a, mu, true, false);
				hi = 1.0;
			}
			
			double u = Rf_runif(lo, hi);
			u = std::max(u, 6.7e-16);
			u = std::min(u, 1.0-6.7e-16);
	
			out = Rcpp::stats::qnorm_1(u, mu, true, false);
			//if inverse cdf fails numerically fall back to slice sampler
			if (!arma::is_finite(out) || out<a || out>b) {
				out = rtnorm_slice(10, mu, 1.0, a, b);
			}
    } else {
    	out=a;
    }
    return(out);
}


double rinvgauss(double mu, double lam) {
	
	double x, y;
	
	y = Rf_rnorm(0.0, 1.0);
	y = y*y;
	x = mu + mu*(mu*y - sqrt(4*mu*lam*y + mu*mu*y*y))/(2*lam);
	if (Rf_runif(0.0, 1.0) > mu/(mu+x)) {
		x = mu*mu/x;
	}
	
	return(x);

}
