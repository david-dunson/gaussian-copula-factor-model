library(roxygen2)
dir='~/Dropbox/R_packages/bfa'
setwd(dir)

#f = paste('~/Documents/svn/bfa/pkg/man/',list.files('~/Documents/svn/bfa/pkg/man/'), sep='')
#lapply(f, unlink)

roxygenize('pkg', 'pkg', copy.package=FALSE, overwrite=TRUE, 
           unlink.target=TRUE)
