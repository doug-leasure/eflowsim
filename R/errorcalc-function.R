errorcalc <- function(s, m, alpha=0.05){
  results <- data.frame(sim=1:s$nsim,
                        k.mean=rep(0,s$nsim), k.base=rep(0,s$nsim), k.high=rep(0,s$nsim), k.var=rep(0,s$nsim), 
                        r.mean=rep(0,s$nsim), r.base=rep(0,s$nsim), r.high=rep(0,s$nsim), r.var=rep(0,s$nsim))
  
  for (i in 1:s$nsim){
    
    # Type I error
    if(s$b1r == 0 & s$b1k == 0 & m[[i]][1,'Pr(>|t|)'] <= alpha) results[i, 'k.mean'] <- results[i, 'r.mean'] <- 1
    if(s$b2r == 0 & s$b2k == 0 & m[[i]][2,'Pr(>|t|)'] <= alpha) results[i, 'k.base'] <- results[i, 'r.base'] <- 1
    if(s$b3r == 0 & s$b3k == 0 & m[[i]][3,'Pr(>|t|)'] <= alpha) results[i, 'k.high'] <- results[i, 'r.high'] <- 1
    if(s$b4r == 0 & s$b4k == 0 & m[[i]][4,'Pr(>|t|)'] <= alpha) results[i, 'k.var'] <- results[i, 'r.var'] <- 1
    
    # Type II error
    if(s$b1r != 0 & m[[i]][1,'Pr(>|t|)'] > alpha) results[i, 'r.mean'] <- 2
    if(s$b2r != 0 & m[[i]][2,'Pr(>|t|)'] > alpha) results[i, 'r.base'] <- 2
    if(s$b3r != 0 & m[[i]][3,'Pr(>|t|)'] > alpha) results[i, 'r.high'] <- 2
    if(s$b4r != 0 & m[[i]][4,'Pr(>|t|)'] > alpha) results[i, 'r.var'] <- 2
    
    if(s$b1k != 0 & m[[i]][1,'Pr(>|t|)'] > alpha) results[i, 'k.mean'] <- 2
    if(s$b2k != 0 & m[[i]][2,'Pr(>|t|)'] > alpha) results[i, 'k.base'] <- 2
    if(s$b3k != 0 & m[[i]][3,'Pr(>|t|)'] > alpha) results[i, 'k.high'] <- 2
    if(s$b4k != 0 & m[[i]][4,'Pr(>|t|)'] > alpha) results[i, 'k.var'] <- 2
    
    # Type III error
    if(s$b1r < 0 & m[[i]][1,'Pr(>|t|)'] <= alpha & m[[i]][1,'Estimate'] > 0) results[i, 'r.mean'] <- 3
    if(s$b2r < 0 & m[[i]][2,'Pr(>|t|)'] <= alpha & m[[i]][2,'Estimate'] > 0) results[i, 'r.base'] <- 3
    if(s$b3r < 0 & m[[i]][3,'Pr(>|t|)'] <= alpha & m[[i]][3,'Estimate'] > 0) results[i, 'r.high'] <- 3
    if(s$b4r < 0 & m[[i]][4,'Pr(>|t|)'] <= alpha & m[[i]][4,'Estimate'] > 0) results[i, 'r.var'] <- 3
    
    if(s$b1r > 0 & m[[i]][1,'Pr(>|t|)'] <= alpha & m[[i]][1,'Estimate'] < 0) results[i, 'r.mean'] <- 3
    if(s$b2r > 0 & m[[i]][2,'Pr(>|t|)'] <= alpha & m[[i]][2,'Estimate'] < 0) results[i, 'r.base'] <- 3
    if(s$b3r > 0 & m[[i]][3,'Pr(>|t|)'] <= alpha & m[[i]][3,'Estimate'] < 0) results[i, 'r.high'] <- 3
    if(s$b4r > 0 & m[[i]][4,'Pr(>|t|)'] <= alpha & m[[i]][4,'Estimate'] < 0) results[i, 'r.var'] <- 3
    
    if(s$b1k < 0 & m[[i]][1,'Pr(>|t|)'] <= alpha & m[[i]][1,'Estimate'] > 0) results[i, 'k.mean'] <- 3
    if(s$b2k < 0 & m[[i]][2,'Pr(>|t|)'] <= alpha & m[[i]][2,'Estimate'] > 0) results[i, 'k.base'] <- 3
    if(s$b3k < 0 & m[[i]][3,'Pr(>|t|)'] <= alpha & m[[i]][3,'Estimate'] > 0) results[i, 'k.high'] <- 3
    if(s$b4k < 0 & m[[i]][4,'Pr(>|t|)'] <= alpha & m[[i]][4,'Estimate'] > 0) results[i, 'k.var'] <- 3
    
    if(s$b1k > 0 & m[[i]][1,'Pr(>|t|)'] <= alpha & m[[i]][1,'Estimate'] < 0) results[i, 'k.mean'] <- 3
    if(s$b2k > 0 & m[[i]][2,'Pr(>|t|)'] <= alpha & m[[i]][2,'Estimate'] < 0) results[i, 'k.base'] <- 3
    if(s$b3k > 0 & m[[i]][3,'Pr(>|t|)'] <= alpha & m[[i]][3,'Estimate'] < 0) results[i, 'k.high'] <- 3
    if(s$b4k > 0 & m[[i]][4,'Pr(>|t|)'] <= alpha & m[[i]][4,'Estimate'] < 0) results[i, 'k.var'] <- 3
    
  }
  return(results)
}
