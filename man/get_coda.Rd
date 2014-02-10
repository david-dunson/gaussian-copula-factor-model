\name{get_coda}
\alias{get_coda}
\title{Get coda object}
\usage{
get_coda(model, loadings = TRUE, scores = FALSE, scale = TRUE,
  positive = FALSE)
}
\arguments{
  \item{model}{a \code{bfa} model}

  \item{loadings}{Return samples of the factor loadings?
  (default is TRUE)}

  \item{scores}{Return samples of the factor scores?
  (default is FALSE)}

  \item{scale}{Return factor loadings on the correlation
  scale? (default is TRUE)}

  \item{positive}{Post-process to enforce positivity
  constraints}
}
\value{
An \code{mcmc} object
}
\description{
Returns an \code{mcmc} object for use in \code{coda}
functions for convergence diagnostics, HPD intervals, etc.
}

