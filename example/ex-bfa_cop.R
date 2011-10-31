require(MASS)
data(UScereal)
UScereal$shelf = factor(UScereal$shelf, ordered=TRUE)
UScereal$vitamins = factor(UScereal$vitamins, ordered=TRUE,
                           levels=c("none", "enriched", "100%"))
fit_cop = bfa_copula(~., data=UScereal[,-1], num.factor=2, nsim=5000, nburn=500, thin=2,
                      normal.dist=rep(0,10), rest=list(c("sugars", 2, "0")),
                      loading.prior="gdp", beta=0.2, keep.scores=T, init.fa=FALSE)
plot_loadings(fit_cop)
biplot(fit_cop, cex=c(0.8, 0.8))

