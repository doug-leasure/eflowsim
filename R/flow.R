#' Create stream flow time series
#' @description Create stream flow time series from historical stream gage data.
#' @param gage character. Gage ID number.
#' @return A data frame.
#' @export

flow <- function(gage='07261500'){
  
  # retrieve data
  d <- dataRetrieval::readNWISdv(gage, '00060')
  d <- d[d$X_00060_00003_cd=='A',c('Date','X_00060_00003')]
  
  # water year
  year <- as.numeric(format(d$Date, format='%Y'))
  month <- as.numeric(format(d$Date, format='%m'))
  year[month>9] <- year[month>9]+1
  d$year <- year
  
  # reformat flow data
  names(d) <- c('date','discharge','year')
  d <- d[,c('year','discharge')]
  
  # identify years with at least 330 days of data
  years <- c()
  for (year in unique(d$year)){
    if(sum(d$year==year) > 330) years <- c(years, year)
  }
  
  if (length(years)==0){
    stop('There were no years with more than 330 days of data.')
  }
  
  # prepare results
  result <- data.frame(row.names=as.character(years), year=years)
  
  for (year in years){
    row <- as.character(year)
    di <- d[d$year==year,'discharge']
    result[row,'mean'] <- mean(di)
    result[row,'base'] <- min(zoo::rollapply(di, width = 7, by = 1, FUN = mean, na.rm=T, align= 'left'))
    result[row,'high'] <- max(zoo::rollapply(di, width = 3, by = 1, FUN = mean, na.rm=T, align= 'left'))
    result[row,'cv'] <- sd(di)/mean(di)  
  }
  row.names(result) <- 1:nrow(result)
  
  if(nrow(result) < 15) warning('This flow record has less than 15 years.')
  
  return(result)  

}

