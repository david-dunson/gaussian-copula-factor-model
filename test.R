set.seed(1)

detach("package:bfa", unload = TRUE)
library.dynam.unload("bfa", system.file(package = "bfa"))
library(bfa)

#Make Datasets
P=7; K=1; N=50
A = matrix(rev(c(0, 0.2, 0.5, 0.8, 1.2 ,2, 3)), nrow=P, ncol=K)
di = 1/sqrt(1 + apply(A^2, 1, sum))
scA = A*di
true = c(scA)

Z = scA %*% matrix(rnorm(K*N), nrow=K) + diag(di)%*%matrix(rnorm(P*N), nrow=P)
d.c.1 = Z

Z = scA %*% matrix(rnorm(K*N), nrow=K) + diag(di)%*%matrix(rnorm(P*N), nrow=P)
d.d.1 = matrix(as.numeric(Z>0), nrow=P)

d.m.1 = d.c.1
d.m.1[c(6,7),] = d.d.1[c(6,7),]

K=3
A = matrix(rnorm(P*K), nrow=P, ncol=K); A[upper.tri(A)] = 0
diag(A) = 2
di = 1/sqrt(1 + apply(A^2, 1, sum))
scA = A*di
true = c(scA)
Z = scA %*% matrix(rnorm(K*N), nrow=K) + diag(di)%*%matrix(rnorm(P*N), nrow=P)
d.c.3 = Z

Z = scA %*% matrix(rnorm(K*N), nrow=K) + diag(di)%*%matrix(rnorm(P*N), nrow=P)
d.d.3 = matrix(as.numeric(Z>0), nrow=P)

d.m.3 = d.c.3
d.m.3[c(6,7),] = d.d.3[c(6,7),]

# Copula factor model, K=1
nsim=2000; nburn=500; thin=1
m.c.1 = bfa(d.c.1, nsim=nsim, nburn=nburn, thin=thin, print.status=1000)

# Gaussian factor model, K=1
nsim=2000; nburn=500; thin=1
m.c.1 = bfa(d.c.1, nsim=nsim, nburn=nburn, thin=thin, print.status=1000,
            normal.dist=rep(1, 7))

# Loadings prior var
m.c.1 = bfa(d.c.1, nsim=nsim, nburn=nburn, thin=thin, print.status=1000,
            loadings.var=2.0)
# Normal
m.c.1.nor = bfa(d.c.1, nsim=nsim, nburn=nburn, thin=thin, print.status=1000,
            loadings.var=2.0, loading.prior='normal')
# GDP
m.c.1.gdp.1 = bfa(d.c.1, nsim=nsim, nburn=nburn, thin=thin, print.status=1000,
            loadings.var=2.0, loading.prior='gdp')
m.c.1.gdp.2 = bfa(d.c.1, nsim=nsim, nburn=nburn, thin=thin, print.status=1000,
            loadings.var=2.0, loading.prior='gdp', alpha=20, beta=0.1)

# Gaussian factor model, K=3
nsim=2000; nburn=500; thin=1
m.c.3 = bfa(d.c.3, num.factor=3, nsim=nsim, nburn=nburn, thin=thin, print.status=1000,
            normal.dist=rep(1, 7), px=FALSE, factor.scales=TRUE, restrict = utri_restrict(7,3))
m.c.3 = bfa(d.c.3, num.factor=3, nsim=nsim, nburn=nburn, thin=thin, print.status=1000,
            normal.dist=rep(1, 7), keep.loadings=FALSE, coda='none', init=TRUE)


detach("package:bfa", unload = TRUE)
library.dynam.unload("bfa", system.file(package = "bfa"))
library(bfa)
rest = utri_restrict(7,3)
rest = rest[sample(1:7,7),]
m.c.3.nor = bfa(d.c.3, num.factor=3, restrict = rest, 
            nsim=nsim, nburn=nburn, thin=thin, print.status=1000,
            loadings.var=2.0, loading.prior='normal')
m.c.3.gdp = bfa(d.c.3, num.factor=3, restrict = rest, 
            nsim=nsim, nburn=nburn, thin=thin, print.status=1000,
            loadings.var=2.0, loading.prior='gdp')