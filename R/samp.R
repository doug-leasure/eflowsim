#' Population samples
#' @description Sample from an eflowsim simulated population time series
#' @param s.dat list. Simulated population. (see ?eflowsim::sim)
#' @param nsites numeric. Number of sample sites.
#' @param sitelength.min numeric. Minimum length of a sample site.
#' @param sitelength.max numeric. Maximum length of a sample site.
#' @param npass.min numeric. Minimum number of removal sampling passes.
#' @param npass.max numeric. Maximum number of removal sampling passes.
#' @param sigmaprop numeric. Standard deviation in proportion of total population extent sampled.
#' @param sigmap numeric. Standard deviation in first-pass detection probabilities.
#' @param b0p numeric. Regression intercept for first-pass detection probabilities.
#' @param b1p numeric. Effect of covariate 1 on first-pass detection probabilities.
#' @param b2p numeric. Effect of covariate 2 on first-pass detection probabilities.
#' @param b3p numeric. Effect of covariate 3 on first-pass detection probabilities.
#' @param b4p numeric. Effect of covariate 4 on first-pass detection probabilities.
#' @param delta numeric. Exponential rate of decline in sampling probabilities in each subsequent pass.
#' @return list.
#' @export

samp <- function(s.dat, nsites=NULL, sitelength.min=100, sitelength.max=100, npass.min=3, npass.max=3, 
                 sigmaprop=1, sigmap=1, b0p=0, b1p=0, b2p=0, b3p=0, b4p=0, delta=0){
  # s.dat=s1; extent=10; nsites=NULL; sitelength.min=100; sitelength.max=100; npass.min=3; npass.max=3; sigmaprop=0.01; sigmap=0.01; b0p=0; b1p=0; b2p=0; b3p=0; b4p=0; delta=0
  
  attach(s.dat, warn.conflicts=F)
  
  if(is.null(nsites)) nsites <- rep(10, nyears)
  
  lag <- max(s.dat$rlag, s.dat$klag)
  if (lag > 1){
    nsites[1:max(0,lag-1)] <- 0
  }
  t1 <- 2 + max(0, lag-1)
  
  ypop <- matrix(nrow=nsim, ncol=nyears)
  Nsite <- ysite <- array(dim=c(nsim, nyears, max(nsites)))
  ypass <- array(dim=c(nsim, nyears, max(nsites), npass.max))
  covs <- array(dim=c(nsim, nyears, max(nsites), 4))
  sitelength <- array(dim=c(nsim, nyears, max(nsites)))
    
  for(s in 1:nsim){
    
    prop <- matrix(nrow=nyears, ncol=max(nsites))
    npasses <- matrix(ncol=max(nsites), nrow=nyears)
    rprop <- matrix(ncol=max(nsites), nrow=nyears)

    Psite <- c()
    PIsite <- pisite <- P <- palpha <- matrix(ncol=max(nsites), nrow=nyears)
    mup <- p <- q <- qprev <- pi <- PI <- array(dim=c(nyears, max(nsites,na.rm=T), npass.max))
    
    qprev[,,1] <- 1
        
    for (t in t1:nyears){
      if(nsites[t] > 0){
        
        # covariates
        covs[s, t, 1:nsites[t], 1:4] <- rnorm(nsites[t]*4, 0, 1)
        
        # site length
        for (j in 1:nsites[t]){
          sitelength[s,t,j] <- runif(1, sitelength.min, sitelength.max)
        }
        
        if(sum(sitelength[s,t,], na.rm=T) > (extent*1000)) {
          sitelength[s,t,] <- sitelength[s,t,]/sum(sitelength[s,t,],na.rm=T)*extent
        }
        
        for (j in 1:nsites[t]){
          prop[t,j] <- sitelength[s,t,j]/(extent*1000)
        }
        
        # number of passes
        for(j in 1:nsites[t]){
          if(npass.min==npass.max) {
            passes <- rep(npass.min, 2)
          } else {
            passes <- seq(npass.min, npass.max)
          }
          npasses[t,j] <- sample(passes,1)
        }
        
        # rprop
        for (j in 1:nsites[t]){
          rprop[t,j] <- boot::inv.logit( rnorm(1, boot::logit(prop[t,j]), sigmaprop) )
        }
        
        if(sum(rprop[t,]) > 1) {
          rprop[t,] <- rprop[t,] / sum(rprop[t,])
          print(paste('warning:  rprop was rescaled for year', t, 'in simulation', s, 'because its sum among sites was greater than 1.'))
        }
        
        
        for (j in 1:nsites[t]){

          ## Regression for detection at pass 1 ##
          mupalpha <- b0p + b1p*covs[s,t,j,1] + b2p*covs[s,t,j,2] + b3p*covs[s,t,j,3] + b4p*covs[s,t,j,4]
          palpha[t,j] <- boot::inv.logit(dnorm(1, mupalpha, sigmap))

          for (m in 1:npasses[t,j]){
            p[t,j,m] <- min(0.99, max(0.01, palpha[t,j] * exp(-delta * (m - 1)) ))

            q[t,j,m] <- 1 - p[t,j,m]
            if(m > 1) { 
              qprev[t,j,m] <- prod(q[t,j,1:(m-1)]) 
            }
            
            pi[t,j,m] <- p[t,j,m] * qprev[t,j,m]
          }
          
          P[t,j] <- 1 - prod(q[t, j, 1:npasses[t,j]])
          
          for (m in 1:npasses[t,j]){ 
            PI[t,j,m] <- pi[t,j,m] / P[t,j]
          }
          
          pisite[t,j] <- rprop[t,j] * P[t,j]
        }
        
        Psite[t] <- sum(pisite[t, 1:nsites[t]])
        for (j in 1:nsites[t]){ 
          PIsite[t,j] <- pisite[t,j] / Psite[t] 
        }
        
        # Generate sampling data
        ypop[s,t] <- rbinom(1, N[s, t], Psite[t])
        
        ysite[s, t, 1:nsites[t]] <- rmultinom(1, ypop[s,t], PIsite[t, 1:nsites[t]])
        
        for (j in 1:nsites[t]){
          Nsite[s,t,j] <- rbinom(1, N[s,t], rprop[t,j])
          ypass[s, t,j,1:npasses[t,j]] <- rmultinom(1, ysite[s,t,j], PI[t,j,1:npasses[t,j]]) 
        }
        
      }#nsites[t]>0
    }#1:nt
  }#1:nsim  
  
  detach(s.dat)
  
  result <- list(Nsite=Nsite, ypop=ypop, ysite=ysite, ypass=ypass, 
                 extent=extent, nsites=nsites, sitelength=sitelength, 
                 npass.min=npass.min, npass.max=npass.max, sigmaprop=1e4, 
                 sigmap=1e4, b0p=0, b1p=0, b2p=0, b3p=0, b4p=0, delta=0.5,
                 covs=covs)

  return(result)
}
