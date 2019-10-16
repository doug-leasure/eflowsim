multipassdat <- function(sampdat){
  nsims <- nrow(sampdat$ypop)
  nyears <- ncol(sampdat$ypop)
  
  covs.out <- data.frame()
  dat.out <- data.frame()
  
  for (sim in 1:nsims){
    for (year in 1:nyears){
      if(sampdat$nsites[year] > 0){
        # covs
        covs.out <- rbind(covs.out, data.frame(site=paste(sim,'.',year,'.',1:sampdat$nsites[year], sep=''),
                                               cov1=sampdat$covs[sim,year,,1], 
                                               cov2=sampdat$covs[sim,year,,2], 
                                               cov3=sampdat$covs[sim,year,,3], 
                                               cov4=sampdat$covs[sim,year,,4]))
        # dat
        for(j in 1:sampdat$nsites[year]){
          nm <- sum(!is.na(sampdat$ypass[sim,year,j,]))
          dat.out <- rbind(dat.out, data.frame(site=rep(paste(sim,'.',year,'.',j, sep=''), nm), pass=1:nm, count=sampdat$ypass[sim,year,j,1:nm]))
        }
      }
    }
  }
  dat.out$site <- as.character(dat.out$site)
  covs.out$site <- as.character(covs.out$site)
  
  return(list(dat=dat.out, covs=covs.out))
}