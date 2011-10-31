#' @include class.R
NULL
#' Plots of factor loadings
#' @param x a bfa model object or a (loadings) matrix
#' @param type either "heatmap" or "credint" (HPD intervals)
#' @param sort_order a permutaion of 1:P (where P is the number of variables) providing sort order
#' for the variables (default is alphabetic)
#' @param color character vector of colors (default Bl-Wh-Rd) (for heatmap)
#' @param scale_name label for the legend (for heatmap)
#' @param prob probability content of HPD interval
#' @return a ggplot2 object
#' @export

plot_loadings <- function(x, type=c("heatmap", "credint"), sort_order=NA, color=NA, scale_name="Loading", 
                          prob=0.95) {
  type = match.arg(type)
  if(type=="heatmap") {
    m = loadings_heat(x, color=color, sort_order=sort_order, scale_name=scale_name)
  } else {
    m = loadings_ci(x, prob=prob, sort_order=sort_order)
  }
  return(m)
}

# Display a heatmap of factor loadings
# @aliases loadings_heat
# @param x a bfa model object or a (loadings) matrix
# @param xlabel x-axis label
# @param ylabel y-axis label
# @param color character vector of colors (default Bl-Wh-Rd)
# @param sort_order a permutaion of 1:P (where P is the number of variables) providing sort order
# for the rows of the loadings matrix
# @param scale_name label for the legend
# @return a ggplot2 object
# @export

loadings_heat <- function(x, color=NA, sort_order=NA, scale_name="Loading") {
  if (class(x)=="bfa") {
    loadings = x$post.loadings.mean
    rownames(loadings) = x$varlabel
  } else {
    loadings = x
  }
  if (is.na(color)) {
    color=c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", 
            "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")
  }
  
	ldf = melt(loadings)
	colnames(ldf) = c("Group", "Factor", "value")
  if (!any(is.na(sort_order))) ldf$Group = factor(ldf$Group, levels=rownames(loadings)[sort_order])
  ldf$Factor = as.factor(ldf$Factor)
  
	breaks = seq(-1,1,by = 0.2)
	cl = colorRampPalette(color)(21)
  lim = c(-1.0, 1)
  p = ggplot(ldf, aes(x=Factor, y=Group, fill=value))
  p = p + scale_fill_gradientn(scale_name, colour=cl, limits = lim, breaks=breaks) 
  p = p + geom_tile()
  
  p = p+ ylab("Variable")+xlab("Factor")
  return(p)

}

# Display HPD intervals of factor loadings
# @param x a bfa model object  
# @param prob probability content of HPD interval
# @export

loadings_ci <- function(x, prob=0.95, sort_order=NA) {

  it=HPDinterval(x, prob=prob)
  it$loadings.mean = data.frame(mean(x)$loadings)
  colnames(it$loadings.mean)=paste("Factor",1:x$K, sep='')
  rownames(it$loadings.mean)=x$varlabel
  d = ldply(it, function(x) data.frame( x, v=rownames(x)) )
  pldf = cast(melt(d), v+variable~...)
  if (any(is.na(sort_order))) {
    sort_order = rev(order(x$varlabel))
  }
  pldf$v = factor(pldf$v, levels=x$varlabel[sort_order])
  m=ggplot(pldf, aes(x=v, y=loadings.mean, ymin=loadings.lower, ymax=loadings.upper))
  m=m+geom_hline(yintercept=0, linetype=2)+xlab("Variable")+ylab("Loading")
  m=m+geom_pointrange()+facet_grid(~variable)+coord_flip()

  return(m)

}

#' Display a biplot
#' @param x A bfa object
#' @param factors Numeric/integer vector giving indices of the factors to plot
#' @param ... Additional arguments to biplot; see \code{?biplot}
#' @method biplot bfa
#' @return Shows a biplot
#' @export
biplot.bfa <- function(x, factors=c(1,2), ...) {
  call_args = list(...)
  if(is.null(call_args$xlabs)) call_args$xlabs = x$obslabel
  if(is.null(call_args$ylabs)) call_args$ylabs = x$varlabel
  call_args$x = t(x$post.scores.mean)[,factors]
  call_args$y = x$post.loadings.mean[,factors]
  do.call(biplot, call_args)
}