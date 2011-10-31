\name{HPDinterval.bfa}
\alias{HPDinterval.bfa}
\title{HPD intervals from a bfa object}
\usage{
  \method{HPDinterval}{bfa} (obj, prob = 0.95, loadings =
  TRUE, scores = FALSE, ...)
}
\arguments{
  \item{obj}{A bfa object}

  \item{prob}{Target probability content (see
  ?HPDinterval)}

  \item{loadings}{Compute interval for the loadings}

  \item{scores}{Compute interval for the scores}

  \item{...}{Ignored}
}
\value{
  A list with elements loadings.lower, loadings.upper,
  scores.lower, scores.upper which are matrices of
  dimension p x k or n x k, or NA's if loadings or scores
  is FALSE
}
\description{
  HPD intervals from a bfa object
}

