\name{bfa_model}
\alias{bfa_model}
\title{Initialize a bfa model}
\usage{
  bfa_model(x, data = NULL, num.factor = 1, restrict = NA,
    normal.dist = NA, center.data = TRUE,
    scale.data = FALSE, init = TRUE, ...)
}
\arguments{
  \item{x}{A formula or matrix}

  \item{data}{The data if x is a formula}

  \item{num.factor}{number of factors}

  \item{restrict}{A matrix or list giving restrictions on
  factor loadings. A matrix should be the same size as the
  loadings matrix. Acceptable values are 0 (identically 0),
  1 (unrestricted), or 2 (strictly positive). List elements
  should be character vectors of the form c('variable',1,
  '>0') where 'variable' is the manifest variable, 1 is the
  factor, and '>0' is the restriction. Acceptable
  restrictions are '>0' or '0'.}

  \item{normal.dist}{A character vector specifying which
  variables should be treated as observed Gaussian}

  \item{center.data}{Center each margin}

  \item{scale.data}{Scale each margin to unit
  mean/variance}

  \item{init}{Initialize the factor loadings}

  \item{...}{ignored}
}
\value{
  An S3 object of class \code{bfa}.
}
\description{
  This function accepts a data matrix \code{D} and
  specified options, returning an S3 object of class sbfac.
}

