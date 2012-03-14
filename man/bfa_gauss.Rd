\name{bfa_gauss}
\alias{bfa_gauss}
\title{Initialize and fit a Gaussian factor model}
\usage{
  bfa_gauss(x, data = NULL, num.factor = 1, restrict = NA,
    center.data = TRUE, scale.data = TRUE, nsim = 0,
    nburn = 0, thin = 1, print.status = 500,
    keep.scores = FALSE, keep.loadings = TRUE,
    loading.prior = c("gdp", "pointmass", "normal"),
    coda = "loadings", ...)
}
\arguments{
  \item{x}{A formula or bfa object.}

  \item{data}{The data if x is a formula}

  \item{num.factor}{Number of factors}

  \item{restrict}{A matrix or list giving restrictions on
  factor loadings. A matrix should be the same size as the
  loadings matrix. Acceptable values are 0 (identically 0),
  1 (unrestricted), or 2 (strictly positive). List elements
  should be character vectors of the form c("variable",1,
  ">0") where 'variable' is the manifest variable, 1 is the
  factor, and ">0" is the restriction. Acceptable
  restrictions are ">0" or "0".}

  \item{center.data}{Center data}

  \item{scale.data}{Scale data}

  \item{nsim}{Number of iterations past burn-in}

  \item{nburn}{Number of initial (burn-in) iterations to
  discard}

  \item{thin}{Keep every thin'th MCMC sample (i.e. save
  nsim/thin samples)}

  \item{print.status}{How often to print status messages to
  console}

  \item{keep.scores}{Save samples of factor scores}

  \item{keep.loadings}{Save samples of factor loadings}

  \item{loading.prior}{Specify point mass ("pointmass",
  default) or normal priors ("normal")}

  \item{coda}{Create \code{mcmc} objects to allow use of
  functions from the \code{coda} package: "all" for
  loadings and scores, "loadings" or "scores" for one or
  the other, or "none" for neither}

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
  Additional parameters: \itemize{ \item loadings.var:
  Factor loading prior variance \item tau.a, tau.b: Gamma
  hyperparameters (scale=1/b) for factor precisions (if
  factor.scales=T) \item rho.a, rho.b: Beta hyperparameters
  for point mass prior \item sigma2.a, sigma2.b: Gamma
  hyperparameters for error precisions \item gdp.alpha,
  gdp.beta: GDP prior parameters }
}

