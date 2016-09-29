#include <float.h>
#include "loadings.h"
#include "bfa_common.h"
//#include <RcppArmadillo.h>
//#include <Rcpp.h>

using namespace arma ;
using namespace Rcpp ;


void sampleLoadings(arma::mat& Z, arma::mat& A, arma::mat& F, 
          Rcpp::NumericVector& sigma2inv,
					Rcpp::NumericVector& error_var_i,
					arma::mat& A_restrict, double A_prior_var, 
					int px, int n, int p, int k, bool positive,
					arma::mat prior_prec) {

  double r, u, v, lo, un, samp, ssq, prior_ssq, scale;
  arma::mat Ap_var, Ap_means, Rf, Apriorv, tt, Ap_var_sm;
  arma::mat U, R;
  arma::vec Ap_mean_sm;
  arma::mat FFt = F*trans(F);
  arma::mat Ftmp, Fsub, disub;
  
  if (prior_prec.n_rows < p) {
    prior_prec = arma::ones<arma::mat>(p,k);
  }
  
  arma::vec sumf2 = arma::sum(arma::pow(F, 2), 1);
  for (int j=0; j<p; j++) {
    arma::mat di = arma::diagmat(prior_prec.row(j));
    //PX step
    if (px>0 && error_var_i(j)<1.0) {
      Ap_var = inv(FFt*sigma2inv(j) + (1/A_prior_var)*di);
    	Ap_means = Ap_var*F*trans(Z.row(j))*sigma2inv(j);
      double ssq = accu(square(trans(Z.row(j)) - trans(F)*Ap_means));
      double prior_ssq = as_scalar(trans(Ap_means) * (1/A_prior_var)*di * Ap_means);
      double scale = 2.0/(ssq + prior_ssq);
      r = Rf_rgamma(0.5*n, scale);
    } else {
      r = 1.0;
    }
    // Sample rows of A
    arma::uvec uu = arma::find(A_restrict.row(j)==0);
    arma::uvec po = arma::find(A_restrict.row(j)==2);
    if (uu.is_empty()) { //No zero restrictions
      Ap_var = inv(FFt*sigma2inv(j) + (1/A_prior_var)*di);
    	Ap_means = Ap_var*F*trans(Z.row(j))*sigma2inv(j);
    	Rf = arma::chol(Ap_var);
			U = my_randn(k, 1);
			arma::vec draw = Rf*U + sqrt(r)*Ap_means;
			
			if (!po.is_empty()) { //Positive restrictions
			  uvec test = find(draw.elem(po)<0);
			  int maxdraw=10000;
			  R_CheckUserInterrupt();
			  while(!test.is_empty() && maxdraw>0) {
			    U = my_randn(k, 1);
			    draw = Rf*U + sqrt(r)*Ap_means;
			    maxdraw--;
			    test = find(draw.elem(po)<0);
			  }
			}
			A.row(j) = trans(draw);
			//for (int i=0; i<k; i++) {
			//	A(j,i) = draw(i);
			//}
    } else { //Zero restrictions
      int lu = uu.n_elem;
      arma::uvec nuu = arma::find(A_restrict.row(j)!=0);
      int nlu = nuu.n_elem;
      Ftmp = F;
      disub = eye(k-lu, k-lu);
      for (int s=0; s<lu; s++) {
        Ftmp.row(uu(s)) = math::inf()*ones(1,n);
      }
      for (int s=0; s<nlu; s++) {
        disub(s,s)= di(nuu(s), nuu(s));
      }
      Fsub = reshape(Ftmp.elem(find(Ftmp<math::inf())), k-lu, n);
      Ap_var_sm = inv(Fsub*trans(Fsub)*sigma2inv(j) + (1/A_prior_var)*disub);
			Ap_mean_sm = Ap_var_sm * Fsub*trans(Z.row(j))*sigma2inv(j);
			U = my_randn(k-lu, 1);
			R = chol(Ap_var_sm);
			vec draw =  R*U + sqrt(r)*Ap_mean_sm;	
			rowvec tmprow = zeros(1,k);
			for (int s=0; s<nlu; s++) {
        tmprow(nuu(s)) = draw(s);
      }
      
      if (!po.is_empty()) { //Positive restrictions
        uvec test = find(tmprow.elem(po)<0);
			  int maxdraw=10000;
			  R_CheckUserInterrupt();
			  while(!test.is_empty() && maxdraw>0) {
			    U = my_randn(k-lu, 1);
			    draw =  R*U + sqrt(r)*Ap_mean_sm;	
			    for (int s=0; s<nlu; s++) {
            tmprow(nuu(s)) = draw(s);
          }
			    maxdraw--;
			    test = find(tmprow.elem(po)<0);
			  }
			}
			A.row(j) = tmprow;
    }
    	//if positive diagonal restriction, sample a_jj | a_j1,...,a_j(j-1)
// 			if (positive) {
// 				tt = Z.row(j) - A.row(j)*F + A(j,j)*F.row(j);
// 				u = arma::accu(F.row(j)%tt)*sigma2inv(j)*sqrt(r);
// 				v = sumf2(j)*sigma2inv(j) + di(j,j)/A_prior_var;
// 				lo = Rcpp::stats::pnorm_1( -u/(sqrt(v)), 0.0, true, false);
// 				un = std::min(Rf_runif(lo, 1.0), 1.0-6.7e-16);
// 				samp = u/v + Rcpp::stats::qnorm_1(un, 0.0, true, false)/sqrt(v);
// 			}
  }
}
