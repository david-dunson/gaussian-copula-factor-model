\name{bfa_mixed}
\alias{bfa_mixed}
\title{Initialize and fit a mixed-scale Gaussian factor model, with probit specifications for the discrete
margins}
\usage{
bfa_mixed(x, data = NULL, num.factor = 1, restrict = NA,
  normal.dist = NA, center.data = TRUE, scale.data = TRUE, nsim = 10000,
  nburn = 1000, thin = 1, print.status = 500, keep.scores = FALSE,
  loading.prior = c("gdp", "pointmass", "normal"), factor.scales = FALSE,
  px = TRUE, coda = "loadings", coda.scale = TRUE, imh.iter = 500,
  imh.burn = 500, ...)
}
\arguments{
  \item{x}{A formula or bfa object.}

  \item{data}{The data if x is a formula}

  \item{num.factor}{Number of factors}

  \item{restrict}{A matrix or list giving restrictions on
  factor loadings. A matrix should be the same size as the
  loadings matrix. Acceptable values are 0 (identically 0),
  1 (unrestricted), or 2 (strictly positive). List elements
  should be character vectors of the form c('variable',1,
  ">0") where 'variable' is the manifest variable, 1 is the
  factor, and ">0" is the restriction. Acceptable
  restrictions are ">0" or "0".}

  \item{normal.dist}{A character vector specifying which
  variables should be treated as observed Gaussian.
  Defaults to all numeric variables if x is a formula.}

  \item{center.data}{Center data (for continuous variables
  - recommended!)}

  \item{scale.data}{Scale data (for continuous variables)}

  \item{nsim}{Number of iterations past burn-in}

  \item{nburn}{Number of initial (burn-in) iterations to
  discard}

  \item{thin}{Keep every thin'th MCMC sample (i.e. save
  nsim/thin samples)}

  \item{print.status}{How often to print status messages to
  console}

  \item{keep.scores}{Save samples of factor scores}

  \item{loading.prior}{Specify GDP ("gdp", default) point
  mass ("pointmass") or normal priors ("normal")}

  \item{factor.scales}{Include a separate scale parameter
  for each factor}

  \item{px}{Use parameter expansion for ordinal variables
  (recommended)}

  \item{coda}{Create \code{mcmc} objects to allow use of
  functions from the \code{coda} package: "all" for
  loadings and scores, "loadings" or "scores" for one or
  the other, or "none" for neither}

  \item{coda.scale}{Put the loadings on the correlation
  scale when creating \code{mcmc} objects}

  \item{imh.iter}{Iterations used to build IMH proposal}

  \item{imh.burn}{Burn-in before collecting samples used to
  build IMH proposal (total burn-in is
  nburn+imh.iter+imh.burn)}

  \item{...}{Prior parameters and other (experimental)
  arguments (see details)}
}
\value{
A \code{bfa} object with posterior samples.
}
\description{
This function performs a specified number of MCMC
iterations and returns an object containing summary
statistics from the MCMC samples as well as the actual
samples if keep.scores or keep.loadings are \code{TRUE}.
Default behavior is to save only the loadings.
}
\details{
Additional parameters: \itemize{ \item loadings.var: Factor
loading prior variance \item tau.a, tau.b: Gamma
hyperparameters (scale=1/b) for factor precisions (if
factor.scales=T). Default is a=b=1 (MV t w/ df=2) \item
rho.a, rho.b: Beta hyperparameters for point mass prior
\item sigma2.a, sigma2.b: Gamma hyperparameters for error
precisions (for numeric variables) \item gdp.alpha,
gdp.beta: GDP prior parameters }
}

