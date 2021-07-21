# cleanup
rm(list=ls()); gc(); cat("\014"); try(dev.off(), silent=T)

# working directory
setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path),'../wd'))

# library
library(multipass)

# load data
readRDS('multipass_data.rds')

# fit multipass model
m1 <- multipass::obsmod(mp1$dat, mp1$covs, n.thin=1, n.burn=5e3, n.samp=5e3)

# traceplots(m1)

q1 <- multipass::quants(m1)

# multipass total site abundance estimates
n <- nrow(mp1$covs)
y <- q[paste('N[',1:n,']',sep=''),'mean']
yup <- q[paste('N[',1:n,']',sep=''),'97.5%']
ylow <- q[paste('N[',1:n,']',sep=''),'2.5%']

# simulated total site abundances
x <- c()
for(i in 1:n){
  ind <- as.numeric(strsplit(as.character(mp$covs$site[i]), "\\.")[[1]])
  x <- c(x, d$Nsite[ind[1], ind[2], ind[3]])
}

plot(y~x, xlab='"real" site abundance', ylab='predicted site abundance')
abline(0,1, col='red')
for(i in 1:length(x)){
  arrows(x0=x[i], x1=x[i], y0=ylow[i], y1=yup[i], length=0)
}




