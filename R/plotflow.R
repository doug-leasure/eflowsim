#' Plot hydrology time-series
#' @description Plot time-series of summary statistics for stream hydrology.
#' @param s list. Simulated population. (see ?eflowsim::sim)
#' @param simid numeric vector.
#' @param type character.
#' @return plot.
#' @export

plotflow <- function(s, simid=NA, type='summary'){
  
  # Setup Plot
  layout(matrix(c(1,2,3,4), nrow=4, ncol=1), heights=c(1.4,1,1,1.4))
  
  # setup data
  x <- 1:s$nyears
  
  if (is.na(simid)) simid <- 1:s$nsim
  
  if(type=='summary'){
    if(length(simid)==1) {
      cov1mean <- cov1up <- cov1low <- s$covs[simid,,1]
      cov2mean <- cov2up <- cov2low <- s$covs[simid,,2]
      cov3mean <- cov3up <- cov3low <- s$covs[simid,,3]
      cov4mean <- cov4up <- cov4low <- s$covs[simid,,4]
      
    } else {
      cov1mean <- apply(s$covs[simid,,1], 2, mean, na.rm=T)
      cov1up <- apply(s$covs[simid,,1], 2, quantile, prob=c(0.975), na.rm=T)
      cov1low <- apply(s$covs[simid,,1], 2, quantile, prob=c(0.025), na.rm=T)
      
      cov2mean <- apply(s$covs[simid,,2], 2, mean, na.rm=T)
      cov2up <- apply(s$covs[simid,,2], 2, quantile, prob=c(0.975), na.rm=T)
      cov2low <- apply(s$covs[simid,,2], 2, quantile, prob=c(0.025), na.rm=T)
      
      cov3mean <- apply(s$covs[simid,,3], 2, mean, na.rm=T)
      cov3up <- apply(s$covs[simid,,3], 2, quantile, prob=c(0.975), na.rm=T)
      cov3low <- apply(s$covs[simid,,3], 2, quantile, prob=c(0.025), na.rm=T)
      
      cov4mean <- apply(s$covs[simid,,4], 2, mean, na.rm=T)
      cov4up <- apply(s$covs[simid,,4], 2, quantile, prob=c(0.975), na.rm=T)
      cov4low <- apply(s$covs[simid,,4], 2, quantile, prob=c(0.025), na.rm=T)
      
    }
  }
  
  # Panel 1
  par(mar=c(1,4.5,4,1))
  if(type=='summary'){
    plot(NA, xlim=range(x), ylim=c(min(cov1low, na.rm=T), max(cov1up, na.rm=T)), 
         xlab=NA, xaxt='n', ylab='Mean Flow',#\nlog-scale', 
         main='Hydrology Time Series',
         cex.lab=1.1, cex.axis=1, cex.main=2)#, log='y')
    lines(x=x, y=cov1mean, lty=1, lwd=1.5)
    lines(x=x, y=cov1up, lty=2, lwd=1)
    lines(x=x, y=cov1low, lty=2, lwd=1)  
  }
  if(type=='spaghetti'){
    plot(NA, xlim=range(x), ylim=c(min(s$covs[,,1], na.rm=T), max(s$covs[,,1], na.rm=T)), 
         xlab=NA, xaxt='n', ylab='Mean Flow',#\nlog-scale', 
         main='Hydrology Time Series',
         cex.lab=1.1, cex.axis=1, cex.main=2)#, log='y')
    
    for(i in simid){
      lines(x=x, y=s$covs[i,,1], col=s$col[i])
    }
    abline(h=0,lty=2)
  }
  axis(side=1, labels=NA)
  
  # Panel 2
  par(mar=c(1,4.5,1,1))  
  if(type=='summary'){
    plot(NA, xlim=range(x), ylim=c(min(cov2low, na.rm=T), max(cov2up, na.rm=T)), 
         ylab='Base Flow', xlab=NA, 
         xaxt='n', cex.lab=1.1, cex.axis=1, lwd=2)
    lines(x=x, y=cov2mean, type='l', lty=1)
    lines(x=x, y=cov2up, type='l', lty=2)
    lines(x=x, y=cov2low, type='l', lty=2)
  }
  if(type=='spaghetti'){
    plot(NA, xlim=range(x), ylim=c(min(s$covs[,,2], na.rm=T), max(s$covs[,,2], na.rm=T)), 
         ylab='Base flow', xlab=NA, 
         xaxt='n', cex.lab=1.1, cex.axis=1, lwd=2)

    for(i in simid){
      lines(x=x, y=s$covs[i,,2], col=s$col[i])
    }
    abline(h=0,lty=2)
  }
  axis(side=1, labels=NA)
  
  # Panel 3
  par(mar=c(1,4.5,1,1))  
  if(type=='summary'){
    plot(NA, xlim=range(x), ylim=c(min(cov3low, na.rm=T), max(cov3up, na.rm=T)), 
         ylab='High Flow', xlab=NA, 
         xaxt='n', cex.lab=1.1, cex.axis=1, lwd=2)
    lines(x=x, y=cov3mean, type='l', lty=1)
    lines(x=x, y=cov3up, type='l', lty=2)
    lines(x=x, y=cov3low, type='l', lty=2)
  }
  if(type=='spaghetti'){
    plot(NA, xlim=range(x), ylim=c(min(s$covs[,,3], na.rm=T), max(s$covs[,,3], na.rm=T)), 
         ylab='High flow', xlab=NA, 
         xaxt='n', cex.lab=1.1, cex.axis=1, lwd=2)
    
    for(i in simid){
      lines(x=x, y=s$covs[i,,3], col=s$col[i])
    }
    abline(h=0,lty=2)
  }
  axis(side=1, labels=NA)
  
  
  # Panel 4
  par(mar=c(4.5,4.5,1,1))  
  if(type=='summary'){
    plot(NA, xlim=range(x), ylim=c(min(cov4low,na.rm=T), max(cov4up,na.rm=T)), 
         ylab='Flow Variability', xlab='Year', 
         cex.lab=1.1, cex.axis=1, lwd=2)
    lines(x=x, y=cov4mean, type='l', lty=1)
    lines(x=x, y=cov4up, type='l', lty=2)
    lines(x=x, y=cov4low, type='l', lty=2)
  }
  if(type=='spaghetti'){
    plot(NA, xlim=range(x), ylim=c(min(s$covs[,,4],na.rm=T), max(s$covs[,,4],na.rm=T)), 
         ylab='Flow Variability', xlab='Year', 
         xaxt='n', cex.lab=1.1, cex.axis=1, lwd=2)
    for(i in simid){
      lines(x=x, y=s$covs[i,,4], col=s$col[i])
    }
    abline(h=0,lty=2)
  }
  axis(side=1)
}