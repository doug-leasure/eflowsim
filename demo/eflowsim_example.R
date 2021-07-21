# cleanup
rm(list=ls()); gc(); cat("\014"); try(dev.off(), silent=T)

# working directory
setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path),'../wd'))

# library
library(eflowsim)

## hydrological data

# option 1: provide flow time series (e.g. example data)
data(flashy.dat)
data(stable.dat)

f1 <- flashy.dat

# option 2: get flow time series from historical USGS stream gage data
f1 <- flow('10249300') #10249300 S Twin River  #10343500 Sagehen  #10329500 Martin Creek  #10318500 Humbolt near Elko  #10317500 Humbolt at Devils Gate (LCT)  #10322000 Maggie Creek (LCT)  #10321000 Humbolt near Carlin (LCT)  #07057500 North Fork (AR Stable) #07261500 Fouche La Favre (AR Flashy)

## population simulation

# settings
parms <- list(covs=f1, nsim=5, nyears=10)
parms <- list(covs=f1, nsim=5, nyears=25, logcovs=T, densityType='phi', extent=10, sigmaR=0.7, b0phi=-0.0025, b1phi=0, b2phi=0, b3phi=0, b4phi=0, b0r=0.4, b1r=0, b2r=0, b3r=0.1, b4r=0, rlag=c(1,1,1,1), dlag=c(0,0,0,0))
parms <- list(covs=f1, nsim=5, nyears=25, logcovs=T, densityType='k', extent=10, sigmaR=0.7, b0k=150, b1k=0, b2k=25, b3k=0, b4k=0, b0r=0.4, b1r=0, b2r=0, b3r=0.1, b4r=0, rlag=c(1,1,1,1), dlag=c(0,0,0,0))

# simulation
s1 <- do.call(sim, parms)

s1$extinct

## plot results

# hydrolgical effects on demographic rates
ploteffects(s1)

# population time-series
plotsim(s1, simid=NA, type='summary')
plotsim(s1, simid=NA, type='spaghetti')

# flow time-series
plotflow(s1, simid=NA, type='summary')
plotflow(s1, simid=NA, type='spaghetti')

## generate population sampling data
d1 <- samp(s1, nsites=rep(c(10),length=s1$nyears), 
           sigmaprop=0.01, sigmap=0.01, 
           b0p=boot::logit(0.4), b1p=0, b2p=0, b3p=0, b4p=0, delta=0.4)

# convert to multipass format:  dat & covs
mp1 <- multipassdat(d1)

# save
saveRDS(mp1, 'multipass_data.rds')

