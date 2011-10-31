\name{reg_samp}
\alias{reg_samp}
\title{Compute samples of regression coefficients}
\usage{
  reg_samp(model, resp.var)
}
\arguments{
  \item{model}{\code{bfa} model object}

  \item{resp.var}{A character vector giving the name(s) of
  the response}
}
\value{
  An array of dimension length(index) x p-length(index) x
  (no. of mcmc samples) with posterior samples of
  regression coefficients
}
\description{
  Compute samples of regression coefficients
}

