rm(list=ls())
gc()

library(eflowsim)

# Set working directory
setwd('c:/RESEARCH/2016 eFlow Simulations/wd')

# Create flow time series from real gage
f1 <- flow('10322000') #10249300 S Twin River  #10343500 Sagehen  #10329500 Martin Creek  #10318500 Humbolt near Elko  #10317500 Humbolt at Devils Gate (LCT)  #10322000 Maggie Creek (LCT)  #10321000 Humbolt near Carlin (LCT)  #07057500 North Fork (AR Stable) #07261500 Fouche La Favre (AR Flashy)

# Create simulated population time-series
# parms <- list(covs=f1, nsim=5, nyears=25, logcovs=T, densityType='phi', extent=10, sigmaR=0.7, b0phi=-0.0025, b1phi=0, b2phi=0, b3phi=0, b4phi=0, b0r=0.4, b1r=0, b2r=0, b3r=0.1, b4r=0, rlag=c(1,1,1,1), dlag=c(0,0,0,0))
parms <- list(covs=f1, logcovs=T, nsim=1, nyears=30, extent=15, sigmaR=0.1, densityType='k', b0k=150, b1k=0, b2k=50, b3k=0, b4k=-25, b0r=0.4, b1r=0, b2r=0, b3r=0.3, b4r=0, rlag=c(1,1,1,1), dlag=c(0,0,0,0))

s1 <- do.call(sim, parms)

s1$extinct

# Plots for simulation
ploteffects(s1)

plotflow(s1, simid=NA, type='summary')
plot(s1, simid=NA, type='summary')

# simids <- c(1, sample(seq(1,s1$nsim), 1))
# simids <- c(2)
# plotflow(s1, simid=simids, type='spaghetti')
# plot(s1, simid=simids, type='spaghetti')

# Generate sample data from each simulation
d1 <- samp(s1, nsites=rep(c(20),length=s1$nyears), 
           sitelength.min=50, sitelength.max=50, npass.min=3, npass.max=3, 
           sigmaprop=0.1, sigmap=0.1, b0p=logit(0.5), b1p=0, b2p=0, b3p=0, b4p=0, delta=0)

# Catch per unit effort
cpue.toggle <- T

if(cpue.toggle){
  cpue <- d1$ysite
  cpue[] <- NA
  
  for (s in 1:nrow(d1$ypop)){
    for (t in 2:ncol(cpue)){
      nsites <- sum(!is.na(d1$ysite[s,t,]))
      for (j in 1:nsites){
        cpue[s, t, j] <- d1$ysite[s, t, j] / (sum(!is.na(d1$ypass[s, t, j,])) * d1$sitelength[s, t, j])
      }
    }
  }
}

# Linear models
m1 <- list()
lag <- 0
y.all <- x1.all <- x2.all <- x3.all <- x4.all <- c()

for (s in 1:nrow(d1$ypop)){
  
  if(cpue.toggle) { y <- apply(cpue[s,-1,], 1, mean)
  } else { y <- apply(d1$Nsite[s,-1,], 1, sum)  }
  
  x1 <- s1$covs[s, -1, 1]
  x2 <- s1$covs[s, -1, 2]
  x3 <- s1$covs[s, -1, 3]
  x4 <- s1$covs[s, -1, 4]
  
  if (lag > 0){
    y <- y[-c(1:lag)]
    x1 <- x1[-c((length(x1)-lag+1):length(x1))]
    x2 <- x2[-c((length(x1)-lag+1):length(x1))]
    x3 <- x3[-c((length(x1)-lag+1):length(x1))]
    x4 <- x4[-c((length(x1)-lag+1):length(x1))]
  }
  
  # m1[[s]] <- summary(lm(y~x1+x2+x3+x4))$coefficients[-1, c(1,4)]
  
  m1[[s]] <- summary(lm(y~x1))$coefficients[-1, c(1,4)]
  m1[[s]] <- rbind( m1[[s]], summary(lm(y~x2))$coefficients[-1, c(1,4)] )
  m1[[s]] <- rbind( m1[[s]], summary(lm(y~x3))$coefficients[-1, c(1,4)] )
  m1[[s]] <- rbind( m1[[s]], summary(lm(y~x4))$coefficients[-1, c(1,4)] )
  
  y.all <- c(y.all, y)
  x1.all <- c(x1.all, x1)
  x2.all <- c(x2.all, x2)
  x3.all <- c(x3.all, x3)
  x4.all <- c(x4.all, x4)
}

par(mfrow=c(2,2))
if (cpue.toggle) {ylab <- 'Catch per effort'
} else { ylab <- 'Total site abundances' }

# lm1 <- lm(y.all ~ x1.all + x2.all + x3.all + x4.all)

plot(y.all~x1.all, xlab='Mean Flow', ylab=ylab, col=rgb(0,0,0,0.2), cex.lab=1.5)
lm1 <- lm(y.all ~ x1.all)
abline(a=lm1$coefficients['(Intercept)'], b=lm1$coefficients['x1.all'], col='red', lwd=2)
pvalue <- round(summary(lm1)$coefficients['x1.all','Pr(>|t|)'], 3)
mtext(paste('p =', pvalue), side=1, adj=0.5, line=-2, col='red', cex=1.5)

plot(y.all~x2.all, xlab='Base Flow', ylab=ylab, col=rgb(0,0,0,0.2), cex.lab=1.5)
lm1 <- lm(y.all ~ x2.all)
abline(a=lm1$coefficients['(Intercept)'], b=lm1$coefficients['x2.all'], col='red', lwd=2)
pvalue <- round(summary(lm1)$coefficients['x2.all','Pr(>|t|)'], 3)
mtext(paste('p =', pvalue), side=1, adj=0.5, line=-2, col='red', cex=1.5)

plot(y.all~x3.all, xlab='High Flow', ylab=ylab, col=rgb(0,0,0,0.2), cex.lab=1.5)
lm1 <- lm(y.all ~ x3.all)
abline(a=lm1$coefficients['(Intercept)'], b=lm1$coefficients['x3.all'], col='red', lwd=2)
pvalue <- round(summary(lm1)$coefficients['x3.all','Pr(>|t|)'], 3)
mtext(paste('p =', pvalue), side=1, adj=0.5, line=-2, col='red', cex=1.5)

plot(y.all~x4.all, xlab='Flow Variability', ylab=ylab, col=rgb(0,0,0,0.2), cex.lab=1.5)  
lm1 <- lm(y.all ~ x4.all)
abline(a=lm1$coefficients['(Intercept)'], b=lm1$coefficients['x4.all'], col='red', lwd=2)
pvalue <- round(summary(lm1)$coefficients['x4.all','Pr(>|t|)'], 3)
mtext(paste('p =', pvalue), side=1, adj=0.5, line=-2, col='red', cex=1.5)

par(mfrow=c(1,1))

#### Assess Error Rates ####
alpha <- 0.05
results <- errorcalc(s1, m1, alpha)

mymean <- function(x, value) { mean(x==value) }

error <- rbind(data.frame(k.mean=s1$b1k, k.base=s1$b2k, k.high=s1$b3k, k.var=s1$b4k, r.mean=s1$b1r, r.base=s1$b2r, r.high=s1$b3r, r.var=s1$b4r),
                         apply(results[,-1], 2, mymean, 1),
                         apply(results[,-1], 2, mymean, 2),
                         apply(results[,-1], 2, mymean, 3))

row.names(error) <- c('Effect', 'Type 1', 'Type 2', 'Type 3')

for(name in names(error)){
  if(error['Effect', name] == 0) error[c('Type 2','Type 3'), name] <- NA
  if(error['Effect', name] > 0) error[c('Type 1'), name] <- NA
}

error <- rbind(error, 1-apply(error[-1,], 2, sum, na.rm=T))
row.names(error)[nrow(error)] <- 'Correct'

print(error)

x <- error['Type 2', c('k.base','k.var','r.high')]
names(x) <- c('k.base.t2','k.var.t2','r.high.t2')
x[,c('k.base.t3','k.var.t3','r.high.t3')] <- error['Type 3', c('k.base','k.var','r.high')]
x[,c('k.mean.t1')] <- error['Type 1', 'k.mean']

labels <- c('k.mean.t1','k.base.t2','k.var.t2','r.high.t2','k.base.t3','k.var.t3','r.high.t3')
x <- as.vector(as.matrix(x[,labels]))
labels <- c('mean','k ~ base','k ~ var','r ~ high','k ~ base','k ~ var','r ~ high')

parmar <- par('mar')
par(mfrow=c(1,1), mar=c(4.5,4.5,1,1))

barplot(x, names.arg=labels, ylim=c(0,1.2), yaxt='n', xlab='Flow Effect', ylab='Error Rate', cex.lab=1.2, cex.names=1.1)
axis(2, at=seq(0, 1, 0.2))

segments(x0=0.2, x1=1.2, y0=1.1)
text(x=0.7, y=1.1, labels='Type I', pos=3, cex=1.1)

segments(x0=1.4, x1=4.8, y0=1.1)
text(x=3.1, y=1.1, labels='Type II', pos=3, cex=1.1)

segments(x0=5, x1=8.4, y0=1.1)
text(x=6.7, y=1.1, labels='Type III', pos=3, cex=1.1)

par(mar=parmar)

## Linear models with CHANGE in population

lag <- 0
y.all <- x1.all <- x2.all <- x3.all <- x4.all <- c()

for (s in 1:nrow(d1$ypop)){
  if(cpue.toggle) { y <- apply(cpue[s,-1,], 1, mean)
  } else { y <- apply(d1$Nsite[s,-1,], 1, sum) }
  
  x1 <- s1$covs[s,-1,1]
  x2 <- s1$covs[s,-1,2]
  x3 <- s1$covs[s,-1,3]
  x4 <- s1$covs[s,-1,4]

  y <- y[-1] - y[-length(y)]

  if(lag==0){
    x1 <- x1[-1]
    x2 <- x2[-1]
    x3 <- x3[-1]
    x4 <- x4[-1]
  }

  if(lag==1){
    x1 <- x1[-length(x1)]
    x2 <- x2[-length(x1)]
    x3 <- x3[-length(x1)]
    x4 <- x4[-length(x1)]
  }

  # m1[[s]] <- summary(lm(y~x1+x2+x3+x4))$coefficients[-1, c(1,4)]

  m1[[s]] <- summary(lm(y~x1))$coefficients[-1, c(1,4)]
  m1[[s]] <- rbind( m1[[s]], summary(lm(y~x2))$coefficients[-1, c(1,4)] )
  m1[[s]] <- rbind( m1[[s]], summary(lm(y~x3))$coefficients[-1, c(1,4)] )
  m1[[s]] <- rbind( m1[[s]], summary(lm(y~x4))$coefficients[-1, c(1,4)] )
  
  y.all <- c(y.all, y)
  x1.all <- c(x1.all, x1)
  x2.all <- c(x2.all, x2)
  x3.all <- c(x3.all, x3)
  x4.all <- c(x4.all, x4)
}

par(mfrow=c(2,2))

if (cpue.toggle) {ylab <- 'Delta catch per effort'
} else { ylab <- 'Delta total site abundances' }

lm1 <- lm(y.all ~ x1.all + x2.all + x3.all + x4.all)

plot(y.all~x1.all, xlab='Mean Flow', ylab=ylab, col=rgb(0,0,0,0.2), cex.lab=1.5)
lm1 <- lm(y.all ~ x1.all)
abline(a=lm1$coefficients['(Intercept)'], b=lm1$coefficients['x1.all'], col='red', lwd=2)
pvalue <- round(summary(lm1)$coefficients['x1.all','Pr(>|t|)'], 3)
mtext(paste('p =', pvalue), side=1, adj=0.5, line=-2, col='red', cex=1.5)

plot(y.all~x2.all, xlab='Base Flow', ylab=ylab, col=rgb(0,0,0,0.2), cex.lab=1.5)
lm1 <- lm(y.all ~ x2.all)
abline(a=lm1$coefficients['(Intercept)'], b=lm1$coefficients['x2.all'], col='red', lwd=2)
pvalue <- round(summary(lm1)$coefficients['x2.all','Pr(>|t|)'], 3)
mtext(paste('p =', pvalue), side=1, adj=0.5, line=-2, col='red', cex=1.5)

plot(y.all~x3.all, xlab='High Flow', ylab=ylab, col=rgb(0,0,0,0.2), cex.lab=1.5)
lm1 <- lm(y.all ~ x3.all)
abline(a=lm1$coefficients['(Intercept)'], b=lm1$coefficients['x3.all'], col='red', lwd=2)
pvalue <- round(summary(lm1)$coefficients['x3.all','Pr(>|t|)'], 3)
mtext(paste('p =', pvalue), side=1, adj=0.5, line=-2, col='red', cex=1.5)

plot(y.all~x4.all, xlab='Flow Variability', ylab=ylab, col=rgb(0,0,0,0.2), cex.lab=1.5)  
lm1 <- lm(y.all ~ x4.all)
abline(a=lm1$coefficients['(Intercept)'], b=lm1$coefficients['x4.all'], col='red', lwd=2)
pvalue <- round(summary(lm1)$coefficients['x4.all','Pr(>|t|)'], 3)
mtext(paste('p =', pvalue), side=1, adj=0.5, line=-2, col='red', cex=1.5)

par(mfrow=c(1,1))

## Error rate ###
results <- errorcalc(s1, m1, alpha)

mymean <- function(x, value) { mean(x==value) }

error <- rbind(data.frame(k.mean=s1$b1k, k.base=s1$b2k, k.high=s1$b3k, k.var=s1$b4k, r.mean=s1$b1r, r.base=s1$b2r, r.high=s1$b3r, r.var=s1$b4r),
               apply(results[,-1], 2, mymean, 1),
               apply(results[,-1], 2, mymean, 2),
               apply(results[,-1], 2, mymean, 3))
row.names(error) <- c('Effect', 'Type 1', 'Type 2', 'Type 3')

for(name in names(error)){
  if(error['Effect', name] == 0) error[c('Type 2','Type 3'), name] <- NA
  if(error['Effect', name] > 0) error[c('Type 1'), name] <- NA
}

error <- rbind(error, 1-apply(error[-1,], 2, sum, na.rm=T))
row.names(error)[nrow(error)] <- 'Correct'

print(error)

x <- error['Type 2', c('k.base','k.var','r.high')]
names(x) <- c('k.base.t2','k.var.t2','r.high.t2')
x[,c('k.base.t3','k.var.t3','r.high.t3')] <- error['Type 3', c('k.base','k.var','r.high')]
x[,c('k.mean.t1')] <- error['Type 1', 'k.mean']

labels <- c('k.mean.t1','k.base.t2','k.var.t2','r.high.t2','k.base.t3','k.var.t3','r.high.t3')
x <- as.vector(as.matrix(x[,labels]))
labels <- c('mean','k ~ base','k ~ var','r ~ high','k ~ base','k ~ var','r ~ high')

parmar <- par('mar')
par(mfrow=c(1,1), mar=c(4.5,4.5,1,1))

barplot(x, names.arg=labels, ylim=c(0,1.2), yaxt='n', xlab='Flow Effect', ylab='Error Rate', cex.lab=1.2, cex.names=1.1)
axis(2, at=seq(0, 1, 0.2))

segments(x0=0.2, x1=1.2, y0=1.1)
text(x=0.7, y=1.1, labels='Type I', pos=3, cex=1.1)

segments(x0=1.4, x1=4.8, y0=1.1)
text(x=3.1, y=1.1, labels='Type II', pos=3, cex=1.1)

segments(x0=5, x1=8.4, y0=1.1)
text(x=6.7, y=1.1, labels='Type III', pos=3, cex=1.1)

par(mar=parmar)