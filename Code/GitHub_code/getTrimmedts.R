require(zoo)

#' Function to remove values in a vector around the median
#' returns a vector shortened by the slice length 
#' 
#' parameters:
#' @param vec:  a univariate vector of NDVI values, 
#'              assumed to have no NA values
#' @param slicelen:  int. The number of items to remove from the vector. 
rmAroundMedian <- function(vec, slicelen){
  # create dataframe with the the values 
  # of the vector and the median of the vector as columns
  med <- data.frame(median = median(unlist(vec)),
                    value = vec)
  
  # add a column 'dist' that calculates 
  # the distance of each value to the median 
  med$dist <- abs(med$median - med$value)
  
  # grab the indices of the values by order of their proximity to the median 
  # and use these to sort the value column 
  removal <- med$value[order(med$dist)] 
  
  # slice the removal column to the desired slice length 
  # and remove these values from the input vector 
  output <- vec[! vec %in% removal[1:slicelen]]
  return(output)
}

#' Function to remove the middle observation in a vector.
#' the function uses the ceiling, 
#' so the middle observation of a vector length == 1 is vec[1]
#' 
#' parameters 
#' @param vec : vector with NDVI values of that year
#' 
#' @param freq: int, the desired frequency to be obtained for DBEST
removeNaByFreq <- function(vec, freq){
  # determine the trim-factor: how many observations should be removed
  trimfact <- length(vec)-freq
  
  # get the indices where the NDVI == NA
  ind <- which(is.na(vec))
  

  #' - If there are no NA-values found in the input data 
  #'   (when the length of the 'ind'-vector is 0),
  #'   it will delete values in the time series that are closest to the median observation 
  #' - If there are NA-values found the number of NA-values 
  #'   necessary to reach the trim-factor are deleted. 
  #'   - If there are not enough NA-values to delete, 
  #'     each additional NA-value is removed around the median observation. 
     
  if(length(ind) > 0){
    if(length(ind) >= trimfact){
      # get a vector with all values minus 
      # the first sequence of redundant Na values 
      vecshort <- vec[-ind[c(1:trimfact)]]
    } else{
      # define how many other values need to be removed 
      trimremainder <-abs(length(ind) - trimfact)
      
      # get a vector with all values minus 
      # the first sequence of redundant Na values 
      vecshort <- unlist(vec[-ind])
      vecshort <- rmAroundMedian(vecshort, trimremainder)
      warning(paste('Not all trimmed values were Na,\nFirst value of the year is:',vecshort[1]))
    }
  }else{
    # throw warning when year has no NA values to remove
    vecshort <- rmAroundMedian(unlist(vec), trimfact)
    warning(paste('No NA values to remove this year.\nFirst value of the year is:',vecshort[1]))
  }
  return(vecshort)
}

#' Function to output yearly vectors at the desired size. 
#' Timeseries of years listed in 'removalyears' 
#' will get trimmed to the desired frequency defined by freq
#' parameters:
#' @param i: 
#'           int. Element of vector to which the function is applied. 
#'           Original vector contains all years (format = "YYYY") 
#'           of which satellite observations are available
#' @param ndvivec: 
#'           univariate vector containing NDVI time series. 
#'           The vector names are assumed to be dates with format "YYYY-MM-DD"
#' @param removalyears:
#'           vector containing years (format: "YYYY") of which the time series should be trimmed
#' @param freq: 
#'           int: desired annual frequency of the time series
trimObsyrs <- function(i, ndvivec, removalyears,freq){
  # slice the ndvi time series to the year of interest
  obs <- ndvivec[which(startsWith(names(ndvivec), prefix = as.character(i)))]
  # evaluate whether the year should be trimmed 
  if(i %in% removalyears){
    return(removeNaByFreq(obs,freq))
  }else{
    return(obs)
  }
}

#' Function to get trimmed and Na-filled time series at a 
#' frequency of 22 observations per year
#' 
#' parameters: 
#' @param i:  index of vector to which function is applied 
#' @param dframe:
#'     data.frame containing NDVI time-series with columns of format "YYYY-MM-DD" 
#' @param lookuptable:
#'     data.frame with "years" and "frequency" containing all years
#'     for which observations are recorded and 
#'     their respective number of observations 
#' @param tsfreq:
#'     int. number specifying the frequency of the output time-series
getTrimmedts <- function(i, dframe, lookuptable, tsfreq, interpolate = TRUE){
  # get ndvi time series from dataframe row
  ndvi <- dframe[i,which(startsWith(colnames(dframe), prefix = "2"))]
  
  if(all(is.na(ndvi))){
    return(as.ts(ndvi))
  }
  # get a vector with the years that have more than 22 observations 
  removalyrs <- lookuptable$Years[which(lookuptable$Frequency > tsfreq)]
  
  # get time series with frequency 22 
  timeseries <- unlist(lapply(lookuptable[,1], FUN = trimObsyrs, ndvivec = ndvi, removalyears = removalyrs, freq = tsfreq ), use.names = TRUE)
  
  # prepare properties for ts object
  thestartyear <- as.numeric(strsplit(names(timeseries[1]),"-")[[1]][1])
  obsnr <- (tsfreq-lookuptable$Frequency[1])+1
  
  # get a ts object
  myts <- ts(data = as.numeric(timeseries), start = c(thestartyear, obsnr), frequency = tsfreq)
  
  # interpolate and extrapolate Na values 
  # na.approx causes collapse of myts into single observation. wur id 8808
  if(interpolate){
    myts <- round(na.approx(myts, rule = 2), digits = 4)
  }
  
  names(myts) <- names(timeseries)
  
  return(myts)
}
