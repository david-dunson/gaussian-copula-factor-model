\name{coef.bfa}
\alias{coef.bfa}
\title{Extract samples of implied regression coefficients from a bfa object}
\usage{
  \method{coef}{bfa} (object, responses,
    scale = attr(object, "type") != "gauss", ...)
}
\arguments{
  \item{object}{A bfa object}

  \item{responses}{A character vector containing one or
  more response variables}

  \item{scale}{Whether to compute regression coefficients
  from factor loadings on the correlation scale;
  recommended if object is a copula or mixed factor model.}

  \item{...}{Ignored}
}
\value{
  An array of dimension length(index) x p-length(index) x
  (no. of mcmc samples) with posterior samples of
  regression coefficients
}
\description{
  Extract samples of implied regression coefficients from a
  bfa object
}

