\name{predict.bfa}
\alias{predict.bfa}
\title{Posterior predictive and univariate conditional posterior predictive distributions}
\usage{
  \method{predict}{bfa} (object, resp.var = NA, cond.vars =
  NA, numeric.as.factor = TRUE, ...)
}
\arguments{
  \item{object}{\code{bfa} model object}

  \item{resp.var}{Either a character vector (length 1) with
  name of the response variable for conditional, or NA for
  draws from the joint posterior predictive.}

  \item{cond.vars}{Conditioning variables; either a list
  like list(X1=val1, X2=val2) with X1, X2 variables in the
  original data frame, or a P length vector with either the
  conditioning value or NA (for marginalized variables).
  Ignored if resp.var is NA}

  \item{numeric.as.factor}{Treat numeric variables as
  ordinal when conditioning}

  \item{...}{Ignored}
}
\value{
  A matrix where each row is either a sample of the
  conditional posterior predictive cdf at each datapoint,
  or a single sample from the joint posterior predictive.
}
\description{
  Posterior predictive and univariate conditional posterior
  predictive distributions, currently implemented only for
  copula models. If resp.var is not NA, returns an estimate
  of the conditional cdf at every observed data point for
  each MCMC iterate. If resp.var is NA, returns draws from
  the joint posterior predictive.
}

