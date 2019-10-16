flow <- function(gage='07261500'){
  d <- readNWISdv(gage, '00060')
  d <- d[d$X_00060_00003_cd=='A',c('Date','X_00060_00003')]
  
  wateryear <- function(date){
    year <- as.numeric(format(date, format='%Y'))
    month <- as.numeric(format(date, format='%m'))
    year[month>9] <- year[month>9]+1
    return(year)
  }
  d$year <- wateryear(d$Date)
  
  names(d) <- c('date','discharge','year')
  d <- d[,c('year','discharge')]
  
  years <- c()
  for (year in unique(d$year)){
    if(sum(d$year==year) > 330) years <- c(years, year)
  }
  
  if (length(years)==0){
    print('There were no years with more than 330 days of data.')
    return()
  } else{
    result <- data.frame(row.names=as.character(years), year=years)
    
    for (year in years){
      row <- as.character(year)
      di <- d[d$year==year,'discharge']
      result[row,'mean'] <- mean(di)
      result[row,'base'] <- min(rollapply(di, width = 7, by = 1, FUN = mean, na.rm=T, align= 'left'))
      result[row,'high'] <- max(rollapply(di, width = 3, by = 1, FUN = mean, na.rm=T, align= 'left'))
      result[row,'cv'] <- sd(di)/mean(di)  
    }
    row.names(result) <- 1:nrow(result)
    
    if(nrow(result) < 15) print('Warning: This flow record has less than 15 years.')
    return(result)  
  }
  
}

