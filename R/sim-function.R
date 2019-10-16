# Generate an eflow population simulation
sim <- function(covs=NULL, logcovs=T, nyears=25, nsim=1, densityType='k', Ninit=NA, growth.max=Inf, pop.max=1e6, extent=14, sigmaR=0.01, b0k=1e4, b1k=0, b2k=1000, b3k=0, b4k=0, b0phi=-1e4, b1phi=0, b2phi=0, b3phi=0, b4phi=0, b0r=1, b1r=0, b2r=0, b3r=0.5, b4r=0, rlag=c(0,0,0,0), dlag=c(0,0,0,0)){
  #nyears=25; covs=NULL; nsim=1; growth.max=Inf; pop.max=1e6; sigmaR=0.01; b0k=1e4; b1k=0; b2k=1e3; b3k=0; b4k=0; b0r=1; b1r=0; b2r=0; b3r=0.5; b4r=-0; lag=c(0,0,0,0)
  
  lag <- max(rlag, dlag)
  
  nyears <- nyears + max(0, lag-1)
  
  if(is.null(covs)) {
    data(flashy.dat)
    covs <- flashy.dat
  }
  
  if(logcovs) {
    covs[,c('mean','base','high','cv')] <- log(covs[,c('mean','base','high','cv')] + 1)
  }
  
  covsim <- array(dim=c(nsim, nyears, 4))
  
  ## Simulate population dynamics
  N <- r <- matrix(nrow=nsim, ncol=nyears)
  
  t1 <- 1 + max(0, lag-1)
  r[,t1] <- b0r
  
  if(densityType == 'k'){
    k <- N
    k[,t1] <- b0k
    if(is.na(Ninit)) {
      N[,t1] <- b0k * extent
    } else {
      N[,t1] <- Ninit
    }
  } else if(densityType == 'phi'){
    phi <- N
    phi[,t1] <- b0phi
    if(is.na(Ninit)) {
      N[,t1] <- max(c(10, b0r/-b0phi), na.rm=T)
    } else {
      N[,t1] <- Ninit
    }
  }
  N[,t1] <- round(N[,t1], 0)
  
  
  for (s in 1:nsim){
    
    if(s == 1){
      covsim[s,,1] <- scale(rep(covs[,'mean'], length=nyears))
      covsim[s,,2] <- scale(rep(covs[,'base'], length=nyears))
      covsim[s,,3] <- scale(rep(covs[,'high'], length=nyears))
      covsim[s,,4] <- scale(rep(covs[,'cv'], length=nyears))
    }else{
      # Sample from hydro data
      rows <- as.character(sample(row.names(covs), size=nyears, replace=T))
      covsim[s,,1] <- scale(covs[rows,'mean'])
      covsim[s,,2] <- scale(covs[rows,'base'])
      covsim[s,,3] <- scale(covs[rows,'high'])
      covsim[s,,4] <- scale(covs[rows,'cv'])
    }
    
    for (t in (t1+1):nyears){
      
      ## Regression for r ##
      r[s, t] <- b0r + b1r*covsim[s,t-rlag[1],1] + b2r*covsim[s,t-rlag[2],2] + b3r*covsim[s,t-rlag[3],3] + b4r*covsim[s,t-rlag[4],4]
      
      ## Regression for density dependence ##
      if(densityType == 'k'){
        r[s,t] <- max(0, r[s, t])
        k[s, t] <- max(1, b0k + b1k*covsim[s,t-dlag[1],1] + b2k*covsim[s,t-dlag[2],2] + b3k*covsim[s,t-dlag[3],3] + b4k*covsim[s,t-dlag[4],4])
        R <- r[s, t] * (1 - (N[s, t-1] / extent) / k[s, t] )
      } 
      if(densityType == 'phi'){
        phi[s, t] <- min(0, b0phi + b1phi*covsim[s,t-dlag[1],1] + b2phi*covsim[s,t-dlag[2],2] + b3phi*covsim[s,t-dlag[3],3] + b4phi*covsim[s,t-dlag[4],4])
        R <- r[s, t] + phi[s,t] * (N[s, t-1] / extent)
      }
      if(R > growth.max){
        R <- growth.max
        print('warning:  R exceeded growth.max')
      }
      
      lambda <- exp( rnorm(1, R, sigmaR) )
      muN <- N[s, t-1] * lambda
      
      if(muN > pop.max){
        muN <- pop.max
        print('warning:  muN exceeded pop.max')
      }
      
      Nfem <- rbinom(1, round(muN), 0.5)
      if(Nfem == muN | Nfem == 0) muN <- 0
      
      N[s, t] <- rpois(1, muN)
    }
  } # End 1:nsim loop
  
  extinct <- mean(N[,nyears] == 0)
  
  if(densityType=='k'){
    result <- list(N=N, r=r, k=k, covs=covsim, extinct=extinct, nsim=nsim, nyears=nyears, 
                   growth.max=growth.max, pop.max=pop.max, densityType=densityType,
                   extent=extent,
                   b0k=b0k, b1k=b1k, b2k=b2k, b3k=b3k, b4k=b4k,
                   b0r=b0r, b1r=b1r, b2r=b2r, b3r=b3r, b4r=b4r, 
                   sigmaR=sigmaR, rlag=rlag, dlag=dlag, col=rainbow(nsim))  
  } else if(densityType == 'phi'){
    result <- list(N=N, r=r, phi=phi, covs=covsim, extinct=extinct, nsim=nsim, nyears=nyears, 
                   growth.max=growth.max, pop.max=pop.max, densityType=densityType,
                   extent=extent,
                   b0phi=b0phi, b1phi=b1phi, b2phi=b2phi, b3phi=b3phi, b4phi=b4phi, 
                   b0r=b0r, b1r=b1r, b2r=b2r, b3r=b3r, b4r=b4r, 
                   sigmaR=sigmaR, rlag=rlag, dlag=dlag, col=rainbow(nsim))
  }
  
  class(result) <- "eflow_sim"
  return(result)
}
