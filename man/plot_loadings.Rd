\name{plot_loadings}
\alias{plot_loadings}
\title{Plots of factor loadings}
\usage{
  plot_loadings(x, type = c("heatmap", "credint"),
    sort_order = NA, color = NA, scale_name = "Loading",
    prob = 0.95)
}
\arguments{
  \item{x}{a bfa model object or a (loadings) matrix}

  \item{type}{either "heatmap" or "credint" (HPD
  intervals)}

  \item{sort_order}{a permutaion of 1:P (where P is the
  number of variables) providing sort order for the
  variables (default is alphabetic)}

  \item{color}{character vector of colors (default
  Bl-Wh-Rd) (for heatmap)}

  \item{scale_name}{label for the legend (for heatmap)}

  \item{prob}{probability content of HPD interval}
}
\value{
  a ggplot2 object
}
\description{
  Plots of factor loadings
}

