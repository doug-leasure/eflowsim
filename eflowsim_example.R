# # Install the package
# install.packages(pkgs='c:/RESEARCH/2016 eFlow Simulations/eflowsim_0.1.0.zip', repos=NULL)
# library(help=eflowsim)

# Load the library
library(eflowsim)

# Set working directory
setwd('c:/RESEARCH/2016 eFlow Simulations/wd')

# Load example data
data(flashy.dat)
data(stable.dat)

# Create flow time series from real gage
f1 <- flow('10249300') #10249300 S Twin River  #10343500 Sagehen  #10329500 Martin Creek  #10318500 Humbolt near Elko  #10317500 Humbolt at Devils Gate (LCT)  #10322000 Maggie Creek (LCT)  #10321000 Humbolt near Carlin (LCT)  #07057500 North Fork (AR Stable) #07261500 Fouche La Favre (AR Flashy)

# Create simulated population time-series
parms <- list(covs=f1, nsim=5, nyears=10)
parms <- list(covs=f1, nsim=5, nyears=25, logcovs=T, densityType='phi', extent=10, sigmaR=0.7, b0phi=-0.0025, b1phi=0, b2phi=0, b3phi=0, b4phi=0, b0r=0.4, b1r=0, b2r=0, b3r=0.1, b4r=0, rlag=c(1,1,1,1), dlag=c(0,0,0,0))
parms <- list(covs=f1, nsim=5, nyears=25, logcovs=T, densityType='k', extent=10, sigmaR=0.7, b0k=150, b1k=0, b2k=25, b3k=0, b4k=0, b0r=0.4, b1r=0, b2r=0, b3r=0.1, b4r=0, rlag=c(1,1,1,1), dlag=c(0,0,0,0))

s1 <- do.call(sim, parms)

s1$extinct

# Plots for simulation
ploteffects(s1)

plot(s1, simid=NA, type='summary')
plot(s1, simid=NA, type='spaghetti')

plotflow(s1, simid=NA, type='summary')
plotflow(s1, simid=NA, type='spaghetti')

# Generate sample data from each simulation
d1 <- samp(s1, nsites=rep(c(10),length=s1$nyears), sigmaprop=0.01, sigmap=0.01, b0p=logit(0.4), b1p=0, b2p=0, b3p=0, b4p=0, delta=0.4)

# Write multipass data:  dat & covs
mp1 <- multipassdat(d1)


##################
#library(multipass)

m1 <- obsmod(mp1$dat, mp1$covs, n.thin=1, n.burn=5e3, n.samp=5e3)

#traceplots(m1)

q1 <- quants(m1)

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




