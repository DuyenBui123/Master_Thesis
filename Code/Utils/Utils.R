# Modified probaV::smoothLoess to plot the result
# This function is by Johannes Eberenz under Expat license:
# https://github.com/johanez/probaV/blob/master/LICENSE.md
smoothLoessPlot <- function (tsx, QC_good = NULL, dates = NULL, 
                             threshold = c(-50, Inf), 
                             res_type = c("distance", "sd_distance", "all", "filled",
                                          "omit", "QC"), 
                             span=0.3, # was 0.3
                             family="gaussian", 
                             threshstat="none", plot=TRUE, ...){
  if (is.null(QC_good)) {
    QC_good <- as.numeric(!is.na(tsx))
  } else {
    QC_good <- as.numeric(QC_good)
  }
  x <- as.numeric(tsx)
  x[QC_good == 0] <- NA
  if(all(is.na(x))){
    warning("Input is all NA")
    return(x)
  }
  if (plot){
    datesforplot <- as.Date(dates, origin = "1970-01-01")
    plot(datesforplot, x, type="o", col = "red",  xlab="", ylab="", ...)
    title(xlab = "Time \U2192", ylab="Reflectance \U2192", line =2)}
  if (is.null(dates)){
    dates <- index(tsx)
  }
  dates <- as.numeric(dates)
  loe <- try(loess(formula = x ~ dates, na.action = "na.omit", span=span, family=family), silent = T)
  if (class(loe) == "try-error"){
    return(x)}
  loe_pred <- try(predict(loe, dates), silent = T)
  if(class(loe_pred) == "try-error"){
    return(x)
  }
  if (plot){
    lines(datesforplot, loe_pred, col="blue")}
  distance <- (loe_pred - x)
  predmae <- mean(abs(distance), na.rm=TRUE)
  predrmse <- sqrt(mean(distance^2, na.rm=TRUE))
  xsd <- sd(x, na.rm=TRUE)
  xmad <- mad(x, na.rm=TRUE)
  if (plot){
    par(cex = 0.9, mai=c(1,0.6,0.3,0.2), xpd = TRUE)
    title(sub = paste0("MAE: ", round(predmae,2), 
                       " RMSE: ", round(predrmse,2), 
                       " sd: ", round(xsd,2), 
                       " mad: ", round(xmad,2), "\n", 
                       "Threshold: ", threshold[1],", ", threshold[2],
                       " Stat used: ", threshstat), line =4, adj = 0)
    legend(x = "bottomright", 
           inset=c(0,-0.56),
           legend = c("NDVI", "Omitted","Loess fit"),
           bty = "n",
           lty = c(0,0, 1),
           pch = c("o","o",""),
           col = c('black', 'red', 'blue'),
           cex = 0.9)}
  threshstat <- switch(threshstat, none=1, sd=xsd, mad=xmad, mae=predmae, rmse=predrmse)
  threshold <- threshold * threshstat
  if (!is.null(threshold)){
    QC_good[distance < threshold[1] & !is.na(distance)] <- 2
    QC_good[distance > threshold[2] & !is.na(distance)] <- 2
  }
  if (class(tsx) == "zoo") {
    tsx <- zoo(cbind(x = as.numeric(tsx), QC_good, filled = loe_pred),
               index(tsx))
    return(tsx)
  } else {
    x_omit <- x
    x_omit[QC_good != 1] <- NA
    if (plot){
      points(datesforplot, x_omit, type="o", col="black")}
    res <- switch(res_type, 
                  all = data.frame(x = as.numeric(tsx),
                                   QC_good = QC_good, 
                                   filled = loe_pred, 
                                   distance = round(distance)),
                  filled = loe_pred, 
                  omit = x_omit, 
                  QC = QC_good, 
                  distance = distance,
                  sd_distance = (distance/sd(x, na.rm = T)))
    return(res)
  }
}


#' Function to obtain remove cloud for time series satellite data
#'
#' @param source is a sf object and a data.frame of time series with seven layers.
#'                It requires to contain band 2 and centroid y column
#'                The names of the time series columns have to
#'                contain day, month and year with format X????.??.??
#'                
#' @param source_path is a string that contains a directory of the "source" parameter
#' @param source2 can be either a sf object and a data.frame of vegetation index (NDVI)
#'                or a string "none", which means that the second source is the same as the first source
#'                The names of the time series columns have to
#'                contain day, month and year with format X????.??.??
#' @return return either filtered source or source2 
filter_time_series <- function(source, source_path, source2) {
  # read the layers of the input database
  inputlayers <- st_layers(source_path)$name
  if(length(inputlayers) > 1){
    # find the first layer that is not blue
    mylayer<- grep("SR_B2", inputlayers, invert = TRUE, value = TRUE, ignore.case = TRUE)[1]
    # read the time series of this layer
    thetimeseries <- st_read(source_path, layer = mylayer)
  } else {
    # read the time series and assume it is the blue band 
    thetimeseries <- st_read(source_path)
  }
  
  #TileList = GetTileList(LoadGlobalRasterPoints())
  
  # get the columns that include the date labels
  ColDates <- colnames(thetimeseries)[startsWith(colnames(thetimeseries), "X2")]
  # get the dates as nummeric values
  Dates <- ColDates %>% 
    str_replace_all("[[:punct:][:cntrl:]]", "") %>% 
    str_extract("[2][0][0-2][0-9][0-1][0-9][0-3][0-9]") %>%
    as.Date( format = "%Y%m%d") %>%
    as.numeric()
  
  
  # find the blue band in the input file to make a BlueMatrix
  if(length(inputlayers) > 1){
    # match for B2 in the layer names, we assume the data is Landsat 8
    BlueSF <- st_read((source_path), layer ="SR_B2")
    # make a vector with Y-coordinate values 
    BlueY <- BlueSF[,grep("y", colnames(BlueSF), ignore.case = TRUE)] 
    gc(TRUE)
    # Make a matrix with values from the blue band 
    BlueMatrix <- BlueSF %>% 
      as.data.frame() %>% 
      select(all_of(ColDates)) %>% 
      as.matrix()
    
  } else {
    # We assume the input data is only the blue band
    # Make a vector with Y-coordinate values 
    BlueY <- thetimeseries[,grep("y", colnames(thetimeseries), ignore.case = TRUE)] 
    # Make a matrix with values from the blue band 
    BlueMatrix <- thetimeseries %>% 
      as.data.frame() %>% 
      select(all_of(ColDates)) %>% 
      as.matrix()
    
  }
  # Get blue mask
  CloudMask <- pbapply(BlueMatrix, 1, smoothLoessPlot, dates=Dates, res_type="omit", 
                       span=0.2, threshstat = "sd", threshold=c(-2, 2), 
                       family="symmetric", plot=FALSE)
  CloudMask <- t(CloudMask)
  
  # Statistics
  mean(is.na(BlueMatrix)) # 60.2% of all data is NA
  mean(is.na(CloudMask)) # 61.4% of all data is now NA, so we only removed 1%
  print(paste("The previous %NA was: ", 
              BlueMatrix %>% is.na() %>% mean() %>% round(digits = 4),
              "the current %NA is: ",
              CloudMask %>% is.na() %>% mean() %>% round(digits = 4)))
  
  
  #' Stop when not all y-coordinate values of InputSF 
  #' are equal to the y-coordinate values of the BlueMatrix (as defined by BlueY)
  #' If a directory is the input for source2, then remove cloud for NDVI dataset
  if (source2 != "none") { # read the NDVI data
    source2 <- sf::st_read(source2)
    stopifnot(all(source2 %>% select(contains("y")) == BlueY))
    
    #' Make a matrix of all columns containing observations. 
    InputMatrix <- source2 %>% 
      as.data.frame() %>% 
      select(starts_with("X20")) %>% 
      as.matrix()
    # Set all observations in InputMatrix to NA that are also NA in the CloudMask
    InputMatrix[is.na(CloudMask)] <- NA
    
    
    # Replace all observations with the filtered observations from InputMatrix
    source2[,grep("X20", colnames(source2), ignore.case = T)] <- InputMatrix
    return(source2)} else {stopifnot(all(source %>% select(contains("y")) == BlueY)) 
      # Filter cloud for source
      
      #' Make a matrix of all columns containing observations. 
      InputMatrix <- source %>% 
        as.data.frame() %>% 
        select(starts_with("X20")) %>% 
        as.matrix()
      
      
      # Set all observations in InputMatrix to NA that are also NA in the CloudMask
      InputMatrix[is.na(CloudMask)] <- NA
      
      
      # Replace all observations with the filtered observations from InputMatrix
      source[,grep("X20", colnames(source), ignore.case = T)] <- InputMatrix
      return(source)}
}

# This function used to create a sample id dataframe and save it as a csv
#' 
#' @param data  a dataframe that contains sample id column
#' @param outut a directory path including the name of output file

#' @return the sample id dataframe
get_sampleid <-  function(data, output) {
  outputname <- paste0(main_dir, output)
  # select only sample id column and convert to integer
  sample_id <- as.integer(data[, names(data) %in% c("sample_id")])
  # convert to dataframe
  sample_id <- as.data.frame(sample_id)
  colnames(sample_id) <- NULL
  write.csv(sample_id, outputname)
  return(sample_id)
}

# This function used to create a year month date dataframe of the time series dates and save it as a csv
#' 
#' @param data  a dataframe that contains at list all the ts dates column. The format of the column names should be X.xxxx.xx.xx
#' @param outut a directory path including the name of output file

#' @return the ymd dataframe
getymd <- function(data, output) {
  outputname_ymd <- paste0(main_dir, output)
  # drop all columns that not dates
  drop <- c("sample_id","id")
  data <- data[, !names(data) %in% drop]
  # extract date names
  date_ts <- colnames(data)
  # extract year
  year_ts <- c()
  for (index in 1:length(date_ts)) { year_ts <- append(year_ts,substring(date_ts[index], 1, 4))
  }
  year_ts_df <- as.data.frame(year_ts)
  
  # extract month
  month_ts <- c()
  for (index in 1:length(date_ts)) { month_ts <- append(month_ts,substring(date_ts[index], 6, 7))
  }
  month_ts_df <- as.data.frame(month_ts)
  
  # extract day
  day_ts <- c()
  for (index in 1:length(day_ts)) { day_ts <- append(day_ts,substring(date_ts[index], 9, 10))
  }
  day_ts_df <- as.data.frame(day_ts)
  # combine year, month, day
  date_ts_df <- cbind(year_ts_df, month_ts_df, day_ts_df)
  colnames(date_ts_df) <- NULL
  write.csv(date_ts_df, outputname_ymd)
  return(date_ts_df)
}

# This function used to create year fraction for each ts date and save it as a csv
#' 
#' @param data  a dataframe that contains at list all the ts dates column. The format of the column names should be X.xxxx.xx.xx
#' @param outut a directory path including the name of output file

#' @return the year fraction dataframe
getfracyear <- function(data, output) {
  drop <- c("sample_id","id")
  data <- data[, !names(data) %in% drop]
  date_ts <- colnames(data)
  # remove the leter X with an empty space
  date_ts <- sub('X','',date_ts)
  # replace all "." to "-"
  date_ts <-  gsub("\\.", "-", date_ts, perl = TRUE)
  # convert year-month-day to year fraction
  t <- c()
  for (index in 1:length(date_ts)) {
    t <- append(t, sprintf("%.11f",decimal_date(as.POSIXlt(date_ts[index]))))
  }
  options(digits=11)
  t <- as.numeric(t)
  t <- as.data.frame(t)
  colnames(t) <- NULL
  outputname_t <- paste0(main_dir, output)
  write.csv(t, outputname_t)
  return(t)
}


