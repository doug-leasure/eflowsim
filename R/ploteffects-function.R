ploteffects <- function(s){
  layout(matrix(1:8, nrow=4, ncol=2), heights=c(1.25,1,1,1))
  
  covnames <- c('Mean Flow', 'Base Flow', 'High Flow', 'Flow Variability')
  
  # Plot covariate effects on r
  for (i in 1:4){
    if(i==1) { 
      par(mar=c(4.5,4,3,1))
      main='Effects on growth rate'
    } else { 
      par(mar=c(4.5,4,1,1)) 
      main=NA}
    
    x <- seq(min(s$covs[,,i]), max(s$covs[,,i]), length=100)
    y <- s$b0r + s[[paste('b',i,'r',sep='')]] * x
    plot(y~x, type='l', xlab=covnames[i], ylab='r', main=main)
  }
  
  # Plot covariate effects on density dependence
  if(s$densityType == 'k'){
    for (i in 1:4){
      if(i==1) { 
        par(mar=c(4.5,4,3,1))
        main='Effects on carrying capacity'
      } else { 
        par(mar=c(4.5,4,1,1)) 
        main=NA}
      
      x <- seq(min(s$covs[,,i]), max(s$covs[,,i]), length=100)
      y <- s$b0k + s[[paste('b',i,'k',sep='')]] * x
      plot(y~x, type='l', xlab=covnames[i], ylab='K', main=main)
    }  
  }
  if(s$densityType == 'phi'){
    for (i in 1:4){
      if(i==1) { 
        par(mar=c(4.5,4,3,1))
        main='Effects on density-dependence'
      } else { 
        par(mar=c(4.5,4,1,1)) 
        main=NA}
      
      x <- seq(min(s$covs[,,i]), max(s$covs[,,i]), length=100)
      y <- s$b0phi + s[[paste('b',i,'phi',sep='')]] * x
      plot(y~x, type='l', xlab=covnames[i], ylab='phi', main=main)
    }
  }
  
}