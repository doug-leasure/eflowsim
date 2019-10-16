plot.eflow_sim <- function(s, simid=NA, type='summary'){
  
  if(is.na(simid)) {simid <- 1:s$nsim}
  
  # setup data
  x <- 1:s$nyears
  
  if(type=='summary'){
    if(length(simid)==1) {
      
      Nmean <- Nup <- Nlow <- s$N[simid,]
      rmean <- rup <- rlow <- s$r[simid,]
      if(s$densityType == 'k') kmean <- kup <- klow <- s$k[simid,]
      if(s$densityType == 'phi') phimean <- phiup <- philow <- s$phi[simid,]
    
    } else {
      
      Nmean <- apply(s$N[simid,], 2, mean, na.rm=T)
      Nup <- apply(s$N[simid,], 2, quantile, prob=c(0.975), na.rm=T)
      Nlow <- apply(s$N[simid,], 2, quantile, prob=c(0.025), na.rm=T)
      
      rmean <- apply(s$r[simid,], 2, mean, na.rm=T)
      rup <- apply(s$r[simid,], 2, quantile, prob=c(0.975), na.rm=T)
      rlow <- apply(s$r[simid,], 2, quantile, prob=c(0.025), na.rm=T)
      
      if(s$densityType == 'k'){
        kmean <- apply(s$k[simid,], 2, mean, na.rm=T)
        kup <- apply(s$k[simid,], 2, quantile, prob=c(0.975), na.rm=T)
        klow <- apply(s$k[simid,], 2, quantile, prob=c(0.025), na.rm=T)
      }
      
      if(s$densityType == 'phi'){
        phimean <- apply(s$phi[simid,], 2, mean, na.rm=T)
        phiup <- apply(s$phi[simid,], 2, quantile, prob=c(0.975), na.rm=T)
        philow <- apply(s$phi[simid,], 2, quantile, prob=c(0.025), na.rm=T)
      }
    }
  }
  
  # Setup Plot
  layout(matrix(c(1,2,3), nrow=3, ncol=1), heights=c(0.45, 0.25, 0.30))
  
  # Top panel: Population size (N)
  par(mar=c(1,5.5,4,1))
  if(type=='summary'){
    plot(NA, xlim=range(x), ylim=c(min(Nlow, na.rm=T), quantile(Nup, probs=c(0.99), na.rm=T)), 
         xlab=NA, xaxt='n', ylab='Population Size (N)',#\nlog-scale', 
         main='Simulated Time Series',
         cex.lab=1.1, cex.axis=1, cex.main=2)#, log='y')
    lines(x=x, y=Nmean, lty=1, lwd=1.5)
    lines(x=x, y=Nup, lty=2, lwd=1)
    lines(x=x, y=Nlow, lty=2, lwd=1)  
  }
  if(type=='spaghetti'){
    plot(NA, xlim=range(x), ylim=c(min(s$N, na.rm=T), quantile(s$N, probs=c(0.99), na.rm=T)), 
         xlab=NA, xaxt='n', ylab='Population Size (N)',#\nlog-scale', 
         main='Simulated Time Series',
         cex.lab=1.1, cex.axis=1, cex.main=2)#, log='y')
    
    colors <- rainbow(s$nsim)
    for(i in simid){
      lines(x=x, y=s$N[i,], col=s$col[i])
    }
  }
  
  axis(side=1, labels=NA)
  
  # Middle panel: Growth Rate (r)
  par(mar=c(1,5.5,1,1))  
  if(type=='summary'){
    plot(NA, xlim=range(x), ylim=c(min(rlow,na.rm=T), max(rup,na.rm=T)), 
         ylab='Population\nGrowth Rate', xlab=NA, 
         xaxt='n', cex.lab=1.1, cex.axis=1, lwd=2)
    abline(h=0)
    lines(x=x, y=rmean, type='l', col='blue', lty=1)
    lines(x=x, y=rup, type='l', col='blue', lty=2)
    lines(x=x, y=rlow, type='l', col='blue', lty=2)
  }
  if(type=='spaghetti'){
    plot(NA, xlim=range(x), ylim=c(min(s$r,na.rm=T), max(s$r,na.rm=T)), 
         ylab='Population\nGrowth Rate', xlab=NA, 
         xaxt='n', cex.lab=1.1, cex.axis=1, lwd=2)
    abline(h=0)
    for(i in simid){
      lines(x=x, y=s$r[i,], col=s$col[i])
    }
  }
  axis(side=1, labels=NA)
  
  # Density-dependence (K or phi)
  par(mar=c(4.5,5.5,1,1))  
  if(type=='summary'){
    if(s$densityType == 'k'){
      plot(NA, xlim=range(x), ylim=c(min(klow,na.rm=T), max(kup,na.rm=T)), 
           ylab='Carrying Capacity', xlab='Year', 
           cex.lab=1.1, cex.axis=1, lwd=2)
      lines(x=x, y=kmean, type='l', col='red', lty=1)
      lines(x=x, y=kup, type='l', col='red', lty=2)
      lines(x=x, y=klow, type='l', col='red', lty=2)  
    }
    if(s$densityType == 'phi'){
      plot(NA, xlim=range(x), ylim=c(min(philow,na.rm=T), max(phiup,na.rm=T)), 
           ylab='Density-dependence (phi)', xlab='Year (t)', 
           cex.lab=1.1, cex.axis=1, lwd=2)
      lines(x=x, y=phimean, type='l', col='red', lty=1)
      lines(x=x, y=phiup, type='l', col='red', lty=2)
      lines(x=x, y=philow, type='l', col='red', lty=2)
    }
    
  }
  if(type=='spaghetti'){
    if(s$densityType=='k'){
      plot(NA, xlim=range(x), ylim=c(min(s$k,na.rm=T), max(s$k,na.rm=T)), 
           ylab='Carrying Capacity', xlab='Year', 
           xaxt='n', cex.lab=1.1, cex.axis=1, lwd=2)
      abline(h=0)
      for(i in simid){
        lines(x=x, y=s$k[i,], col=s$col[i])
      }  
    }
    if(s$densityType=='phi'){
      plot(NA, xlim=range(x), ylim=c(min(s$phi,na.rm=T), max(s$phi,na.rm=T)), 
           ylab='Density-dependence (phi)', xlab='Year (t)', 
           xaxt='n', cex.lab=1.1, cex.axis=1, lwd=2)
      abline(h=0)
      for(i in simid){
        lines(x=x, y=s$phi[i,], col=s$col[i])
      }
    }
  }
  axis(side=1)
}