# @nord
.updateScores <- function (Z_, A_, F_) 
.Call("updateScoresC", Z_, A_, F_, PACKAGE = "bfa")

# @nord
 .updateRho <- function (rho_, A_, rhoa_, rhob_) 
.Call("updateRho", rho_, A_, rhoa_, rhob_, PACKAGE = "bfa")

# @nord
#.updateSparseLoadingsJ( SEXP Z_, SEXP A_, SEXP F_, SEXP tauinv_, SEXP rho_, SEXP A_restrict_, SEXP pnz_ )
 .updateSparseLoadings <- function (Z_, A_, F_, tauinv_, rho_, A_restrict_, pnz_) 
.Call("updateSparseLoadingsJ", Z_, A_, F_, tauinv_, rho_, A_restrict_, pnz_, PACKAGE = "bfa")

#SEXP updateZ( SEXP Z_, SEXP Ra_, SEXP maxes_, SEXP argsorts_, SEXP A_, SEXP F_ ){

# @nord
 .updateZ <- function (Z_, Ra_, maxes_, argsorts_, A_, F_) {
.Call("updateZ", Z_, Ra_, maxes_, argsorts_, A_, F_,  PACKAGE = "bfa")
return()
}

# @nord
.updateZcut <- function (Z_, Ra_, maxes_, argsorts_, A_, F_) {
.Call("updateZcut", Z_, Ra_, maxes_, argsorts_, A_, F_,  PACKAGE = "bfa")
return()
}

# @nord
#SEXP MCMCstep( SEXP Z_, SEXP A_, SEXP F_, SEXP tauinv_, SEXP rho_ , SEXP Ra_, SEXP maxes_, SEXP argsorts_, SEXP priors_, SEXP nsim_, SEXP nburn_, SEXP thin_)
.MCMCstep <- function (Z_, A_, F_, tauinv_, rho_, Ra_, maxes_, argsorts_, A_restrict_, priors_, nsim_, nburn_, thin_, printstatus_, keep.scores, keep.loadings, more_args)
.Call("MCMCstep", Z_, A_, F_, tauinv_, rho_, Ra_, maxes_, argsorts_, A_restrict_, priors_, nsim_, nburn_, thin_, printstatus_, keep.scores, keep.loadings, more_args, PACKAGE = "bfa")

#' Initialize and fit a Gaussian factor model
#'
#' This function performs a specified number of MCMC iterations and
#' returns an object containing summary statistics from the MCMC samples
#' as well as the actual samples if keep.scores or keep.loadings are \code{TRUE}.
#' Default behavior is to save only the loadings. 
#' 
#' Additional parameters:
#' \itemize{
#' \item loadings.var: Factor loading prior variance
#' \item tau.a, tau.b: Gamma hyperparameters (scale=1/b) for factor precisions (if factor.scales=T)
#' \item rho.a, rho.b: Beta hyperparameters for point mass prior
#' \item sigma2.a, sigma2.b: Gamma hyperparameters for error precisions
#' \item gdp.alpha, gdp.beta: GDP prior parameters
#' }
#' 
#' @param x A formula or bfa object. 
#' @param data The data if x is a formula
#' @param num.factor Number of factors
#' @param restrict A matrix or list giving restrictions on factor loadings. A matrix should be the 
#' same size as the loadings matrix. Acceptable values are 0 (identically 0), 1 (unrestricted), 
#' or 2 (strictly positive). List elements should be character vectors of the form c("variable",1, ">0")
#' where 'variable' is the manifest variable, 1 is the factor, and ">0" is the restriction. Acceptable
#' restrictions are ">0" or "0".
#' @param center.data Center data
#' @param scale.data  Scale data
#' @param nsim Number of iterations past burn-in
#' @param nburn Number of initial (burn-in) iterations to discard
#' @param thin Keep every thin'th MCMC sample (i.e. save nsim/thin samples)
#' @param print.status How often to print status messages to console
#' @param keep.scores Save samples of factor scores
#' @param keep.loadings Save samples of factor loadings
#' @param loading.prior Specify point mass ("pointmass", default) or normal priors ("normal") 
#' @param coda Create \code{mcmc} objects to allow use of functions from the 
#' \code{coda} package: "all" for loadings and scores, "loadings" or "scores" for one or the
#' other, or "none" for neither
#' @param ... Prior parameters and other (experimental) arguments (see details)
#' @return A \code{bfa} object with posterior samples.
#' @export

bfa_gauss <- function(x, data=NULL, num.factor=1, restrict=NA, 
                center.data=TRUE, scale.data=TRUE, nsim=0, nburn=0, thin=1,
                print.status=500, keep.scores=FALSE, keep.loadings=TRUE,
                loading.prior=c("gdp", "pointmass", "normal"), 
                coda="loadings", ...) {
  loading.prior = match.arg(loading.prior)
  fr = model.frame(x, data=data, na.action=NULL)
  d = dim(fr)
  normal.dist = rep(1, d[2])
  m = .bfa(x, data=data, num.factor=num.factor, restrict=restrict, normal.dist=normal.dist, 
       center.data=center.data, scale.data=scale.data, nsim=nsim, nburn=nburn, thin=thin,
       print.status=print.status, keep.scores=keep.scores, keep.loadings=keep.loadings,
       loading.prior=loading.prior, factor.scales=TRUE, px=FALSE, coda=coda, 
       imh=FALSE, ...) 
  attr(m, "type") <- "gauss"
  return(m)
                
}

#' Initialize and fit a copula factor model
#'
#' This function performs a specified number of MCMC iterations and
#' returns an object containing summary statistics from the MCMC samples
#' as well as the actual samples if keep.scores or keep.loadings are \code{TRUE}.
#' Default behavior is to save only the loadings. 
#' 
#' Additional parameters:
#' \itemize{
#' \item loadings.var: Factor loading prior variance
#' \item tau.a, tau.b: Gamma hyperparameters (scale=1/b) for factor precisions (if factor.scales=T)
#' \item rho.a, rho.b: Beta hyperparameters for point mass prior
#' \item gdp.alpha, gdp.beta: GDP prior parameters
#' }
#' 
#' @param x A formula or bfa object. 
#' @param data The data if x is a formula
#' @param num.factor Number of factors
#' @param restrict A matrix or list giving restrictions on factor loadings. A matrix should be the 
#' same size as the loadings matrix. Acceptable values are 0 (identically 0), 1 (unrestricted), 
#' or 2 (strictly positive). List elements should be character vectors of the form c('variable',1, ">0")
#' where 'variable' is the manifest variable, 1 is the factor, and ">0" is the restriction. Acceptable
#' restrictions are ">0" or "0".
#' @param normal.dist A character vector specifying which variables should be treated as observed 
#' Gaussian. Defaults to none (a completely semiparametric copula model)
#' @param center.data Center continuous variables
#' @param scale.data  Scale continuous variables
#' @param nsim Number of iterations past burn-in
#' @param nburn Number of initial (burn-in) iterations to discard
#' @param thin Keep every thin'th MCMC sample (i.e. save nsim/thin samples)
#' @param print.status How often to print status messages to console
#' @param keep.scores Save samples of factor scores
#' @param keep.loadings Save samples of factor loadings
#' @param loading.prior Specify point mass ("pointmass", default) or normal priors ("normal") 
#' @param factor.scales Include a separate scale parameter for each factor
#' @param px Use parameter expansion for ordinal/copula/mixed factor models (strongly recommended)
#' @param coda Create \code{mcmc} objects to allow use of functions from the 
#' \code{coda} package: "all" for loadings and scores, "loadings" or "scores" for one or the
#' other, or "none" for neither
#' @param coda.scale Put the loadings on the correlation scale when creating \code{mcmc} objects
#' @param imh Use Independence Metropolis-Hastings step for discrete margins. If FALSE, use the 
#' semiparametric extended rank likelihood
#' @param imh.iter Iterations used to build IMH proposal
#' @param imh.burn Burn-in before collecting samples used to build IMH proposal (total burn-in is nburn+imh.iter+imh.burn)
#' @param ... Prior parameters and other (experimental) arguments (see details)
#' @return A \code{bfa} object with posterior samples.
#' @export
#' @examples \dontrun{
#' require(MASS)
#' data(UScereal)
#' UScereal$shelf = factor(UScereal$shelf, ordered=TRUE)
#' UScereal$vitamins = factor(UScereal$vitamins, ordered=TRUE,
#'                            levels=c("none", "enriched", "100%"))
#' obj = bfa_copula(~., data=UScereal[,-1], num.factor=2, nsim=10000, nburn=1000, thin=4,
#'                      rest=list(c("sugars", 2, "0")), loading.prior="gdp", keep.scores=T, 
#'                      print.status=2500)
#' plot_loadings(obj)
#' #plot_loadings returns ggplot objects you can tweak
#' m = plot_loadings(obj, type='credint')
#' print(m)
#' print(m+opts(title="HPD intervals (p=0.95)")+theme_bw())
#' biplot(obj, cex=c(0.8, 0.8))
#' plot(get_coda(obj)) 
#' plot(get_coda(obj, loadings=F, scores=T))
#' hpd.int = HPDinterval(obj, scores=T)
#' 
#' #sample from posterior predictive
#' ps = predict(obj)
#' 
#' m=ggplot(UScereal, aes(x=calories, y=sugars))+geom_point(position='jitter', alpha=0.5)
#' m=m+stat_density2d(data=ps, aes(x=calories, y=sugars, color = ..level..), geom='contour')
#' print(m)
#' 
#' m=ggplot(UScereal, aes(x=calories))+geom_histogram()
#' m=m+stat_density(data=ps, aes(x=calories, y=..count..), color='red',fill=NA, adjust=1.3)
#' m=m+facet_grid(shelf~.)
#' print(m)
#' 
#' #we can compute conditional dist'n estimates directly as well by supplying cond.vars
#' cond.vars=list(shelf=1)
#' out = predict(obj, resp.var="calories", cond.vars=cond.vars)
#' plot(sort(unique(UScereal$calories)), apply(out, 2,mean), type='s')
#' lines(sort(unique(UScereal$calories)), apply(out, 2, quantile, 0.05), type='s', lty=2)
#' lines(sort(unique(UScereal$calories)), apply(out, 2,quantile, 0.95), type='s', lt=2)
#' lines(ecdf(UScereal$calories[UScereal$shelf==1]), col='blue')
#' text(400, 0.1, paste("n =", sum(UScereal$shelf==1)))
#' 
#' out2 = predict(obj, resp.var="calories", cond.vars=list(shelf=2))
#' out3 = predict(obj, resp.var="calories", cond.vars=list(shelf=3))
#' plot(sort(unique(UScereal$calories)), apply(out, 2,mean), type='s')
#' lines(sort(unique(UScereal$calories)), apply(out2, 2,mean), type='s', lty=2)
#' lines(sort(unique(UScereal$calories)), apply(out3, 2,mean), type='s', lty=3)

#' }


bfa_copula <- function(x, data=NULL, num.factor=1, restrict=NA, normal.dist=NA, 
                center.data=TRUE, scale.data=TRUE, nsim=0, nburn=0, thin=1,
                print.status=500, keep.scores=FALSE, keep.loadings=TRUE,
                loading.prior=c("gdp", "pointmass", "normal"), 
                factor.scales=FALSE, px=TRUE,
                coda="loadings", coda.scale=TRUE, imh=FALSE, imh.iter=500,
                imh.burn=500, ...) {
  loading.prior = match.arg(loading.prior)
  if (class(x)!='bfa') {
    fr = model.frame(x, data=data, na.action=NULL)
    d = dim(fr)
    if(any(is.na(normal.dist))) {
      normal.dist = rep(0, d[2])
    }
  }
  m = .bfa(x, data=data, num.factor=num.factor, restrict=restrict, normal.dist=normal.dist, 
       center.data=center.data, scale.data=scale.data, nsim=nsim, nburn=nburn, thin=thin,
       print.status=print.status, keep.scores=keep.scores, keep.loadings=keep.loadings,
       loading.prior=loading.prior, factor.scales=factor.scales, px=px, coda=coda, 
       coda.scale=coda.scale, imh=imh, imh.iter=imh.iter, imh.burn=imh.burn, ...)
  attr(m, "type") <- "copula"
  return(m)
}

#' Initialize and fit a mixed-scale Gaussian factor model, with probit specifications for the discrete
#' margins
#'
#' This function performs a specified number of MCMC iterations and
#' returns an object containing summary statistics from the MCMC samples
#' as well as the actual samples if keep.scores or keep.loadings are \code{TRUE}.
#' Default behavior is to save only the loadings. 
#' 
#' Additional parameters:
#' \itemize{
#' \item loadings.var: Factor loading prior variance
#' \item tau.a, tau.b: Gamma hyperparameters (scale=1/b) for factor precisions (if factor.scales=T)
#' \item rho.a, rho.b: Beta hyperparameters for point mass prior
#' \item sigma2.a, sigma2.b: Gamma hyperparameters for error precisions (for numeric variables)
#' \item gdp.alpha, gdp.beta: GDP prior parameters
#' }
#' 
#' @param x A formula or bfa object. 
#' @param data The data if x is a formula
#' @param num.factor Number of factors
#' @param restrict A matrix or list giving restrictions on factor loadings. A matrix should be the 
#' same size as the loadings matrix. Acceptable values are 0 (identically 0), 1 (unrestricted), 
#' or 2 (strictly positive). List elements should be character vectors of the form c('variable',1, ">0")
#' where 'variable' is the manifest variable, 1 is the factor, and ">0" is the restriction. Acceptable
#' restrictions are ">0" or "0".
#' @param normal.dist A character vector specifying which variables should be treated as observed 
#' Gaussian. Defaults to all numeric variables if x is a formula.
#' @param center.data Center data
#' @param scale.data  Scale data
#' @param nsim Number of iterations past burn-in
#' @param nburn Number of initial (burn-in) iterations to discard
#' @param thin Keep every thin'th MCMC sample (i.e. save nsim/thin samples)
#' @param print.status How often to print status messages to console
#' @param keep.scores Save samples of factor scores
#' @param keep.loadings Save samples of factor loadings
#' @param loading.prior Specify point mass ("pointmass", default) or normal priors ("normal") 
#' @param factor.scales Include a separate scale parameter for each factor
#' @param px Use parameter expansion for ordinal variables (recommended)
#' @param coda Create \code{mcmc} objects to allow use of functions from the 
#' \code{coda} package: "all" for loadings and scores, "loadings" or "scores" for one or the
#' other, or "none" for neither
#' @param coda.scale Put the loadings on the correlation scale when creating \code{mcmc} objects
#' @param imh.iter Iterations used to build IMH proposal
#' @param imh.burn Burn-in before collecting samples used to build IMH proposal (total burn-in is nburn+imh.iter+imh.burn)
#' @param ... Prior parameters and other (experimental) arguments (see details)
#' @return A \code{bfa} object with posterior samples.
#' @export

bfa_mixed <- function(x, data=NULL, num.factor=1, restrict=NA, normal.dist=NA, 
                center.data=TRUE, scale.data=TRUE, nsim=0, nburn=0, thin=1,
                print.status=500, keep.scores=FALSE, keep.loadings=TRUE,
                loading.prior=c("gdp", "pointmass", "normal"), 
                factor.scales=FALSE, px=TRUE,
                coda="loadings", coda.scale=TRUE, imh.iter=500,
                imh.burn=500, ...) {
  loading.prior = match.arg(loading.prior)
  m = .bfa(x, data=data, num.factor=num.factor, restrict=restrict, normal.dist=normal.dist, 
       center.data=center.data, scale.data=scale.data, nsim=nsim, nburn=nburn, thin=thin,
       print.status=print.status, keep.scores=keep.scores, keep.loadings=keep.loadings,
       loading.prior=loading.prior, factor.scales=factor.scales, px=px, coda=coda, 
       coda.scale=coda.scale, imh=TRUE, imh.iter=imh.iter, imh.burn=imh.burn, ...)
  attr(m, "type") <- "mixed"
  return(m)
}

.bfa <- function(x, data=NULL, num.factor=1, restrict=NA, normal.dist=NA, 
                center.data=TRUE, scale.data=TRUE, nsim=0, nburn=0, thin=1,
                print.status=500, keep.scores=FALSE, keep.loadings=TRUE,
                loading.prior="pointmass", factor.scales=FALSE, px=TRUE,
                coda="loadings", coda.scale=TRUE, imh=FALSE, imh.iter=500,
                imh.burn=500, ...) {
  if (class(x) != 'bfa') {
    model = bfa_model(x, data, num.factor, restrict, normal.dist, center.data,
                    scale.data, ...)
  } else {
    model=x
  }
  
  if(imh) {
    cat("Building proposal...")
    model = fit_bfa(model, imh.iter, imh.burn, 1, print.status,
                   keep.scores, keep.loadings, loading.prior,
                   factor.scales, px, coda, coda.scale, save_max=TRUE, quiet=TRUE, ...)

    P = model$P
    means <- vector("list", P)
    covs  <- vector("list", P)
    precs <- vector("list", P)
    
    levs = rep(0, P)
    
    for (ind in 1:P) {
      ind.m = which(model$maxes[ind,]==99999.0) - 1
      levs[ind] = ind.m
      if (model$error_var_i[ind]>0) {
        prop.cov = matrix(1)
        prop.prec = matrix(1)
        prop.mean = matrix(1)
      } else if (ind.m > 1) {
        ps = t(model$post.cutpoints[ind,1:ind.m,])
        sample.cov = 2.4*2.4*cov(as.matrix(ps))/ind.m
        prop.cov = chol(sample.cov)
        prop.prec = solve(sample.cov)
        prop.mean = apply(as.matrix(ps), 2, mean)
      } else {
        ps = t(model$post.cutpoints[ind,1:ind.m,])
        
        prop.cov = sqrt(var(c(ps)))
        prop.prec = 1/prop.cov
        prop.mean = mean(ps)
      }
      covs[[ind]] = prop.cov
      precs[[ind]] = prop.prec
      means[[ind]] = prop.mean
    }
    cat("Done.\n")
    model = fit_bfa(model, nsim, nburn, thin, print.status,
                   keep.scores, keep.loadings, loading.prior,
                   factor.scales, px, coda, coda.scale, imh=1, 
                   prop.cov=covs, prop.prec=precs, prop.mean=means, 
                   nlevels=matrix(levs), df=100,  ...)
  } else {
    model = fit_bfa(model, nsim, nburn, thin, print.status,
                   keep.scores, keep.loadings, loading.prior,
                   factor.scales, px, coda, coda.scale,  ...)
  }
  return(model)

}


# Perform MCMC model fitting for a bfa model
#
# This function performs a specified number of MCMC iterations and
# returns an object containing summary statistics from the MCMC samples
# as well as the actual samples if keep.scores or keep.loadings are \code{TRUE}.
# Default behavior is to save only the loadings. 
# 
# Prior parameters:
# loadings.var: Factor loading prior variance
# tau.a, tau.b: Gamma hyperparameters (scale=1/b) for factor precisions (if factor.scales=T)
# rho.a, rho.b: Beta hyperparameters for point mass prior
# sigma2.a, sigma2.b: Gamma hyperparameters for error precisions
# gdp.alpha, gdp.beta: GDP prior parameters
# 
# @param model an object of type bfa, as returned by bfa(data)
# @param nsim number of iterations past burn-in
# @param nburn number of initial (burn-in) iterations to discard
# @param thin keep every thin'th MCMC sample (i.e. save nsim/thin samples)
# @param print.status how often to print status messages to console
# @param keep.scores save samples of factor scores
# @param keep.loadings save samples of factor loadings
# @param loading.prior Specify point mass ("pointmass", default) or normal priors ("normal") 
# @param factor.scales Include a separate scale parameter for each factor
# @param px Use parameter expansion for ordinal/copula/mixed factor models (recommended)
# @param coda create \code{mcmc} objects to allow use of functions from the 
# \code{coda} package: "all" for loadings and scores, "loadings" or "scores" for one or the
# other, or "none" for neither
# @param coda.scale put the loadings on the correlation scale when creating \code{mcmc} objects
# @param ... Prior parameters and other (experimental) arguments (see details)
# @return The S3 \code{bfa} object \code{model}, now with posterior samples/summaries.

fit_bfa <- function(model, nsim, nburn, thin=1, print.status=500,
                   keep.scores=FALSE, keep.loadings=TRUE, loading.prior="pointmass",
                   factor.scales=FALSE, px=TRUE, coda="loadings", coda.scale=TRUE,  ...) {
  
  more_args = list(...)
  if (is.null(more_args$init)) more_args$init=TRUE
  more_args$px = as.numeric(px)

  rhoa  = ifelse(is.null(more_args$rho.a),     1.0, more_args$rho.a)
  rhob  = ifelse(is.null(more_args$rho.b),     1.0, more_args$rho.b)
  s     = ifelse(is.null(more_args$s),         1.0, more_args$s)
  taua  = ifelse(is.null(more_args$tau.a),     1.0, more_args$tau.a)
  taub  = ifelse(is.null(more_args$tau.b),     1.0, more_args$tau.b)
  s2a   = ifelse(is.null(more_args$sigma2.a),  0.5, more_args$sigma2.a)
  s2b   = ifelse(is.null(more_args$sigma2.b),  1.0, more_args$sigma2.b)
  alpha = ifelse(is.null(more_args$gdp.alpha), 3.0, more_args$gdp.alpha)
  beta  = ifelse(is.null(more_args$gdp.beta),  1.0, more_args$gdp.beta)
  bp    = ifelse(is.null(more_args$bp),       10.0, more_args$gdp.alpha)
  bq    = ifelse(is.null(more_args$bq),        5.0, more_args$gdp.beta)
  

  model$priors = list(rhoa=rhoa, rhob=rhob, s=s, taua=taua, taub=taub, 
                      s2a=s2a, s2b=s2b, alpha=alpha, beta=beta,
                      bp=bp, bq=bq)
  
  if (model$nsim == 0 && more_args$init==TRUE && more_args$px>0) {
    model = .fit(model, max(1000, nburn/10), 0, thin=1, factor.scales=factor.scales,
                 print.status=max(1000, nburn/2)+1,
                 keep.scores=FALSE, keep.loadings=FALSE, 
                 loading.prior=loading.prior, coda="none", 
                 save_max=0, px=0, more_args$error_var_i, quiet=TRUE )
  }
  
  model = .fit(model, nsim, nburn, thin=thin, print.status=print.status,
                   keep.scores=keep.scores, keep.loadings=keep.loadings, 
                   loading.prior=loading.prior, coda=coda, coda.scale=coda.scale,
                   px=more_args$px, ...)
  
}

########################################################
# Actual model fitting fcn
########################################################

.fit <- function(model, nsim, nburn, thin=1, print.status=500,
                   keep.scores=FALSE, keep.loadings=TRUE, loading.prior="pointmass",
                   coda="loadings", coda.scale=TRUE,  ...) {
      
  K = model$K; P = model$P; N = model$N
  
  more_args = list(...)
  
  more_args$loadings_var  = ifelse(is.null(more_args$loadings.var), 1.0, more_args$loadings.var)
  more_args$factor_scales = ifelse(is.null(more_args$factor.scales), FALSE, more_args$factor.scales)
  
  if (loading.prior=="normal")     more_args$method = 0
  if (loading.prior=="pointmass")  more_args$method = 1
  if (loading.prior=="gdp")        more_args$method = 2
  if (loading.prior=="beta")       more_args$method = 3
  
  if(is.null(more_args$imh)) more_args$imh = 0
  if(is.null(more_args$px)) more_args$px = 1
  if(is.null(more_args$df)) more_args$df = 100
  if(is.null(more_args$save_max)) more_args$save_max = 0
  if(is.null(more_args$positive)) more_args$positive = FALSE
  
  more_args$save_max = as.numeric(more_args$save_max)
  more_args$error_var_i = model$error_var_i

	model$nsim  = nsim
	model$nburn = nburn
	model$thin  = thin

	sim = .MCMCstep(model$ldata, model$loadings, model$scores, model$tauinv, 
		              model$rho, model$ranks, model$maxes, model$argsorts, model$loadings.restrict,
		              model$priors, nsim, nburn, thin, print.status, keep.scores, keep.loadings, more_args)
  
  if (keep.scores)   dim(sim$Fp)=c(K,N,nsim/thin)
  if (keep.loadings) dim(sim$Ap)=c(P,K,nsim/thin)
  if (!is.null(more_args$save_max) && more_args$save_max>0) dim(sim$maxp) = c(dim(model$maxes), nsim/thin)	
 
  dim(sim$Fp.mean) = c(K,N)
  dim(sim$Fp.var)  = c(K,N)
  dim(sim$Ap.mean) = c(P,K)
  dim(sim$Ap.var)  = c(P,K)
  dim(sim$pnz)     = c(P,K)

  model$post.loadings.mean = sim$Ap.mean
  model$post.loadings.var  = sim$Ap.var
  model$post.loadings      = sim$Ap
	
  model$post.scores.mean   = sim$Fp.mean
  model$post.scores.var    = sim$Fp.var
  model$post.scores        = sim$Fp 
  
  model$post.sigma2        = 1/sim$sigma2inv
  model$post.loadings.prob = 1.0-sim$pnz
  
 if (!is.null(more_args$save_max) && more_args$save_max>0) {
    model$post.cutpoints = sim$maxp
 }
 
  if (coda!="none") {
    if (coda=="all") {
      model$loadings.mcmc = get_coda(model, loadings=TRUE, scores=FALSE, scale=coda.scale)
      model$scores.mcmc   = get_coda(model, loadings=FALSE, scores=TRUE)
    } else if (coda=="loadings") {
      model$loadings.mcmc = get_coda(model, loadings=TRUE, scores=FALSE, scale=coda.scale)
    } else if (coda=="scores") {
      model$scores.mcmc   = get_coda(model, loadings=FALSE, scores=TRUE)
    }
  }

  return(model)
}


# Imputation on the original scale of the data
#
# This function imputes data on the original scale given a current value of the latent
# continuous data, using the empirical cdf.
#
# @param D Observed data
# @param Z Latent data (standardized)
# @return A vector of length \code{length(D[is.na(D)])} containing the
# imputed values
# @export


#impute <- function(D, Z) {
#  out = NA
#  for (i in 1:dim(D)[1]) {
#    im = quantile(D[i,], pnorm(Z[i,is.na(Z[i,]), type=1, na.rm=TRUE)
#    out = c(out, im)
#  }
#}



