# load packages
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(here, bfast, strucchangeRcpp)

#source utils 
source(here("GitHub_code", "Utils.R"))
library(strucchangeRcpp)
library(bfast)
source(here("GitHub_code/cglops-change-detection/src/utils/enable_fast_bfast.r"))
source(here("GitHub_code/cglops-change-detection/src/bfast-cal/plotting.r"))


#' Functions in this script are adapted from '02-detectbreaks.r' by Dainius Masiliunas, 
#' at: https://github.com/GreatEmerald/cglops-change-detection/blob/master/src/bfast-cal/02-detectbreaks.r#L32

# Calculate dates from number of elements in the input
# Returns a ts object
GetTS <-  function(data, frequency=22, start=2013, years=8)
{
  stopifnot(is.vector(data)) # Univariate only
  # 8-daily frequency is 46, 16-daily is 23, we have 10 years of data
  if (is.null(frequency))
    frequency <- length(data)/years
  return(ts(data, start=start, frequency = frequency))
}


#' BFAST Monitor break detection function
#' @param InputTS matrix or ts, input time series
#' @param monitor_years vector of starts of monitoring periods
#' @param monitor_length Length of monitoring period, in years
#' @param cloud_threshold Minimal number of observations to run the algorithm
#' @param TargetYears Time of breaks as reported by reference, for plotting.
#' @param quiet Suppress output
#' @param plot Plot diagnostic plots, one per monitor_years
#' @param ... Additional arguments for bfastmonitor
#' 
#' Function outputs a data.frame with columns 
#' "Breakpoint" containing the date of a detected break
#' and "Magnitude" containing the magnitude of a detected break. 
TestBFASTMonitor <- function(InputTS, monitor_years=2016:2018, monitor_length=1.25,
                             cloud_threshold=42, TargetYears=NULL, quiet=FALSE, plot=FALSE, ...)
{
  # The input should be a single row of a matrix.
  if (!is.ts(InputTS))
    InputTS <- GetTS(InputTS) # Convert into a ts object
  
  Observations <- sum(!is.na(window(InputTS, end=min(monitor_years)+monitor_length)))
  if (!quiet)
    print(paste("Observations for point:", Observations))
  
  Result <- data.frame("Breakpoint" = NULL, "Magnitude" = NULL)
  if (Observations < cloud_threshold)
  {
    if (!quiet)
      print("too cloudy")
    Result <- data.frame("Breakpoint" = "cloudy", "Magnitude" = "cloudy")
    return(Result)
  }
  
  for (StartYear in monitor_years)
  {
    # Cut time series to length
    ShortTS <- window(InputTS, end=StartYear+monitor_length)
    BM <- try(bfastmonitor(ShortTS, StartYear, ...))
    
    if ("try-error" %in% class(BM)) {
      print("An error has occurred:")
      print(BM)
      next # Assume no break
    }
    
    if (plot)
    {
      plot(BM, xlab="Time", ylab="NIRv * 255")
      abline(v=TargetYears, col="blue")
    }
    # It can happen that the break is detected early next year, then we just discard those
    if (!is.na(BM$breakpoint) && BM$breakpoint < StartYear+1)
      Result <- rbind(Result, c(BM$breakpoint, BM$magnitude))
  }
  if (length(Result) == 0){
    Result <- data.frame("Breakpoint" = -9999, "Magnitude" = 0)
  }
  colnames(Result) <- c("Breakpoint", "Magnitude")
  return(Result)
}

#' Run BFAST Lite over the time series and get the detected breaks.
#'
#' @param scsig significance value at which the structure change test is considered successful
#' @param scrange range over which the sctest should be run
#' @param sctype type of sctest (?efp)
#' @param maginterval interval (% time series) over which to compute breakpoint magnitudes
#' @param magcomponent components (possibly multiple for which to calculate the magnitudes (trend/harmonsin1/harmoncos3/...)
#' @param magstat statistic for magnitude thresholding (diff/RMSD/MAD/MD)
#' @param magthreshold threshold above which the breaks are kept. Breaks lower than the threshold are discarded
#' @param coefcomponent component (one!) for coefficient difference calculation
#' @param coefthresholds min and max, if coefficient difference is within the range, the breakpoint will be discarded.
#' @param plot Whether to call plot.bfast0n() on the output
#' @param quiet Suppress print output
#' @param order Harmonic order
#' @param formula Formula passed to sctest() and bfastpp()
#' @param TargetYears Year fractions of expected breaks for plotting
#' @param seasonfreq Multiplier for how many bins to use in a season. 1 means one bin per observation, 0.5 means two observations per bin.
#'
#' @return Fractional years of all detected breaks, or FALSE if none detected,
#' or NA if not enough observations/error in running the function.
UseBFASTLite <-  function(InputTS, scrange=NULL, scsig=0.05, breaks="LWZ",
                          sctype="OLS-MOSUM", maginterval=0.1, magcomponent="trend",
                          magstat="RMSD", magthreshold=-Inf, coefcomponent="trend",
                          coefthresholds=c(0, 0), plot=FALSE, quiet=FALSE, order=3,
                          formula=response ~ harmon + trend, TargetYears=NULL,
                          seasonfreq=0.5, breaknumthreshold=Inf, altformula=NULL, ...)
{
  # The input should be a single row of a matrix.
  if (!is.ts(InputTS)){
    InputTS <- GetTS(InputTS)} # Convert into a ts object
  
  # define the total number of valid (non-NA) observations
  Observations <- sum(!is.na(InputTS))
  
  if (!quiet){
    print(paste("Observations for point:", Observations))}
  
  # Set h to number of observations per year, i.e. frequency of time series
  h <- ceiling(frequency(InputTS)) 
  
  # if the total number of observations is higher than two times the number of observations per year,
  # we start pre-processing, otherwise we say the time-series contains too much clouds. 
  if (Observations > h*2) {
    bpp <- bfastpp(InputTS, order=order, sbins=seasonfreq) # Preprocess the ts into a data.frame
    
    # Fix season when it's not an integer
    myseason <- as.numeric(as.character(bpp$season)) # Deparse season again
    
    if (!all(is.na(myseason))){
      # If we failed to deparse, then we're using a fixed bfastpp already, no need to do anything
      bpp$season <- cut(myseason, frequency(InputTS)*seasonfreq, ordered_result = TRUE) # Rebin all
    } 
    
    # # Run the sctest first to determine whether to run bfastlite
    # if (!is.null(scrange)) {
    #   SC <- sctest(efp(formula,
    #                    data=bpp[bpp$time > scrange[1] & bpp$time < scrange[2],],
    #                    h=h/Observations, sctype=type))$p.value > scsig
    #   if (!is.null(SC) && !is.na(SC) && SC) 
    #     return(Result <- data.frame("Breakpoint" = NA, "Magnitude" = NA)) # No break detected, return FALSE
    # }
    
    #' Break was detected, so run bfastlite
    #' Run bfastlite using the forumula, the pre-processed regressions curves from `bpp`
    #' and h, the number of observations per year
    bp <- try(breakpoints(formula, data=bpp, h=h))
    
    if ("try-error" %in% class(bp)) {
      print("An error has occurred:")
      print(bp)
      return(Result <- data.frame("Breakpoint" = "Error", "Magnitude" = "Error", "Coefficients" = "Error", "Rsq" = "Error"))
    }
    
    #' now, obtain the timing of the breakpoints, 
    #' using the information criterion of `breaks`, which is `"LWZ"` by default
    bpOptim <- breakpoints(bp, breaks=breaks) # Get breakpoint time
    
    # get the Rsquared 
    rsquared <- bpOptim$r.squared
    
    if (length(bpOptim$breakpoints) > 0 && !all(is.na(bpOptim$breakpoints))) {
      # Are we overfitting? If so, try running with altformula
      if (length(bpOptim$breakpoints) > breaknumthreshold)
      {
        return(MODDetectBreaks(InputTS=InputTS, scrange=scrange, scsig=scsig,
                               breaks=breaks, sctype=sctype, maginterval=maginterval,
                               magcomponent=magcomponent, magstat=magstat, magthreshold=magthreshold,
                               coefcomponent=coefcomponent, coefthresholds = coefthresholds, plot=plot,
                               quiet=quiet, order=order, formula = altformula, TargetYears=TargetYears,
                               seasonfreq=seasonfreq, breaknumthreshold = Inf, ...))
      }
      
      # Get magnitudes and coefficients of each break
      bpMag <- magnitude(bp, breaks=breaks, interval=maginterval, component=magcomponent)$Mag
      coefs <- coef(bp, breaks=breaks)
      bpCoef <- coefs[2:nrow(coefs),coefcomponent] - coefs[1:nrow(coefs)-1,coefcomponent]
      names(bpCoef) <- NULL
    } else {
      bpMag <- NULL
      bpCoef <- NULL
    }
    
    if (plot){
      plot.bfast0n(bp, bpp, breaks=breaks, bpMag=bpMag, ...)
      abline(v=TargetYears, col="red")
    }
    
    if (all(is.na(bpOptim$breakpoints))){ # bfastlite didn't find any breaks, return FALSE
      return(Result <- data.frame("Breakpoint" = -9999, "Magnitude" = 0, "Coefficients" = 0, "Rsq" = rsquared))}
    if (!quiet){
      print("Breakpoints before filtering:")
      print(bpp[bpOptim$breakpoints, "time"])
      print("Magnitudes of each breakpoint:")
      print(bpMag)
      print(paste("Differences in coefficients of", coefcomponent))
      print(bpCoef)
    }
    
    # Convert the output into dates of break
    Result <- data.frame("Breakpoint" = bpp[bpOptim$breakpoints, "time"], 
                         "Magnitude" = bpMag, # Now all output of bpMag will be included in a new column
                         # and appended to the "Magnitude"-label automatically, 
                         # in case you want to output only 
                         # selected magnitude results, use column indexing, like: bpMag[,3]
                         "Coefficients" = bpCoef, "Rsq" = rsquared) 
    
    if (!is.null(bpMag) && !is.null(bpCoef)){
      # Keep breakpoints that are above the threshold of magnitude
      MagFilter <- abs(bpMag[,magstat]) > magthreshold
      # Keep breakpoints where the coefficient difference is big enough
      CoefFilter <- bpCoef < min(coefthresholds) | bpCoef > max(coefthresholds)
      Result <- Result[MagFilter & CoefFilter]
      Result$Rsq <- rsquared
      #' The current code using the indexing of Result,
      #' might snap if more than one break is found.In that case, 
      #' try to select a column in the Magnitude output, like here below: 
      
      #data.frame("Breakpoint" = bpp[bpOptim$breakpoints, "time"][MagFilter & CoefFilter],
      #           "Magnitude" = bpMag[MagFilter][3],
      #           "Coefficients" = bpCoef[CoefFilter])
    }
    if (length(Result) < 1){
      return(Result <- data.frame("Breakpoint" = NA, "Magnitude" = NA, "Coefficients" = NA, "Rsq" = NA))} # If we filtered out all results, return FALSE again
    if (plot){ # For ones that got kept, plot black on top
      abline(v=Result, col="black")}
    return(Result)
    
  } else {
    if (!quiet){
      print("too cloudy")}
    return(Result <- data.frame("Breakpoint" = NA, "Magnitude" = NA, "Coefficients" = NA, "Rsq" = NA))
  }
}

#' Plot the time series and results of bfast0n
#' 
#' @param bp     breakpoints object from strucchange::breakpoints()
#' @param bpp    data.frame output from bfastpp()
#' @param breaks Number of breaks or optimal break selection method, see strucchange::breakpoints()
#' @param bpMag  Break magnitudes, as returned by magnitude(bp)
#' @param ...    Other parameters to pass to plot()
plot.bfast0n = function(bp, bpp, breaks, bpMag=NULL, ...)
{print("im in plot but at the begining")
  # Plot the original time series
  plot(response~time, data=bpp, type="l", ...)
  # Plot the fitted model in green
  lines(fitted(bp, breaks=breaks)~bpp$time, col="green")
  
  # Get the requested breaks
  bpOptim = breakpoints(bp, breaks=breaks)
  #title(sub=paste("RSS:", round(bpOptim$RSS)))
  
  if (length(bpOptim$breakpoints) > 0 && !all(is.na(bpOptim$breakpoints))) {
    bpTimes = bpp[bpOptim$breakpoints, "time"]
    bpY = bpp[bpOptim$breakpoints, "response"]
    print(bpY)
    abline(v=bpTimes, col="blue") # Detected breakpoints in blue
    # If magnitudes requested, plot whiskers denoting magnitude
    arrows(bpTimes, bpY-bpMag[,"RMSD"], bpTimes, bpY+bpMag[,"RMSD"], length=0.05, angle=90, code=3, col="blue")
  }
}
# 3b) t-test
TestMODttest = function(i, ChangedVITS, TargetYears=AllTargetYears, sig=0.05)
{
  Point1TS = ChangedVITS[i, ]
  TargetYear = TargetYears[i]
  if (is.na(TargetYear))
    TargetYear = 2016
  start = as.Date(paste(TargetYear, "01", "01", sep="-"))
  end = start+366
  Observations = sum(!is.na(Point1TS))
  print(paste("Observations for point:", Observations))
  if (Observations > 42) {
    Point1Mean = mean(Point1TS[dates < start], na.rm=TRUE)
    P1Remainder = Point1TS[dates > start & dates < end]
    result = t.test(P1Remainder, mu=Point1Mean)$p.value
    return(result < sig)
    # Might be useful to t-test the other params too
    #P1TS = bfastts(Point1TS, dates, "10-day")
    #P1PP = bfastpp(P1TS, order=3)
    
  } else {
    print("too cloudy")
    return(NA)
  }
}
# This function used to aggregate the output of satellite data after running it on BFAST lite
#' 
#' @param sample_ids     a unique list of samle ids from satellite data
#' @param SR_breakpoints    the output breakpoints of each band of the satellite data. It should at least include:
#'                            sample_id, Breakpoint, and Magnitude columns
#' @param threshold a number that determine whether there is any breakpoints or not                          
#' @return agg_BFASTlite_output a data frame that contain the largest breakpoint magnitude among all bands if there are breakpoints
aggregatefun <- function(sample_ids, SR_breakpoints, threshold) {
  flooryear <- c() # create an empty vector to store the floor year or base year of each breakpoint
  temp_df_gr <- data.frame() # create an empty data frame for temporary storage
  agg_BFASTlite_output <- data.frame() # create an empty data frame to store the output. 
  for (sample_id in sample_ids) { # reset flooryear and temp_data frame for each iteration of a sample id
    flooryear <- c()
    temp_df <- data.frame()
    x_temp_df <- data.frame()
    for (layer in seq_along(SR)) {
      x <- SR_breakpoints[[layer]][SR_breakpoints[[layer]]$sample_id == sample_id,]
      x$flooryear <- "none" # create a new column for the satellite dataset to store floor year of each breakpoint
      flooryear <- append(flooryear, floor(x$Breakpoint)) # append a floor each of each observation in the flooryear vector
    }
    if (all(is.na(flooryear)) == TRUE | all(unique(flooryear) == -9999.000) == TRUE){
      # if all floor years of all bands are either NA or -9999.000, there is no breakpoint. Therefore no need more investigation
      agg_BFASTlite_output <- rbind(agg_BFASTlite_output, x) # append the observation directly to the output data frame
      
    } else if (as.character(names(which.max(table(flooryear)))) == "NA" | as.numeric(names(which.max(table(flooryear)))) == -9999 ) {
      # if the floor years of all the observation contain both NA, -9999.000 and a real numeric number, and NA or -9999.000 are more 
      # dominant then there is no breakpoint occur for that sample id
      for (layer in seq_along(SR_breakpoints)) {
        x <- SR_breakpoints[[layer]][SR_breakpoints[[layer]]$sample_id == sample_id,]
        x$flooryear <- "none"
        if (all(unique(x$Breakpoint) == -9999) == TRUE | all(is.na(x$Breakpoint)) == TRUE) {
          x_temp_df <- rbind(x_temp_df, x)
        }
        
      }
      agg_BFASTlite_output <- rbind(agg_BFASTlite_output, x_temp_df[1,]) 
    } else {
      
      flooryear <- as.data.frame(table(flooryear)) # convert collected floor year of all observation in to a frequency table data frame
      flooryear$flooryear <- as.numeric(as.character(flooryear$flooryear)) # convert factor into numeric
      # extract only floor year that are in the preiod of 2016:2019
      flooryear<- flooryear[flooryear$flooryear %in% flooryear$flooryear[flooryear$flooryear >= 2016 &
                                                                           flooryear$flooryear %in% flooryear$flooryear[flooryear$flooryear <= 2019]],]
      # check whether there is any floor year that occur more than the threshold in order to determine the dominant years that contains threholds
      if(any(flooryear$Freq >= threshold)) { # only floor year that has the frequency larger than the threshold
        selected_flooryear <- flooryear[flooryear$Freq >=threshold,]$flooryear
        for (layer in seq_along(SR_breakpoints)) { # extract all the observations of all bands that contain the dominant floor years
          x <- SR_breakpoints[[layer]][SR_breakpoints[[layer]]$sample_id == sample_id,]
          x$flooryear <- floor(x$Breakpoint)
          x<- x[x$flooryear %in% selected_flooryear,] # selct only observation that has the floor year the same as in the slected floor year of all the layers
          temp_df <- rbind(temp_df,x) # bind the observation into a temporary dataframe
          
        }
        # find the maximum value of magnitude for each group floor year and turn the result into dataframe
        temp_df_gr <- temp_df %>% group_by(flooryear) %>%
          filter(Magnitude == max(Magnitude))%>%
          ungroup() %>%
          as.data.frame()
        
      } else { # there is no floor year that has frequency larger than the threshold, no change will be assigned to that sample id
        obs_nochange <- x[1,]
        obs_nochange$Breakpoint <- -9999.000
        obs_nochange$Magnitude <- 0.0000
        obs_nochange$Coefficients <- 0.0000
        obs_nochange$Rsq <- 0.0000000
        obs_nochange$Magnitude.before <- NA
        obs_nochange$Magnitude.after <- NA
        obs_nochange$Magnitude.diff <- NA
        obs_nochange$Magnitude.RMSD <- NA
        obs_nochange$Magnitude.MAD <- NA
        obs_nochange$Magnitude.MD <- NA
        agg_BFASTlite_output <- rbind(agg_BFASTlite_output, obs_nochange)
      }
      # add the output into the data frame output
      agg_BFASTlite_output <- rbind(agg_BFASTlite_output, temp_df_gr)
    }
    
  }
  return(agg_BFASTlite_output)
}


# Function used for calibrate BFAST lite with different parameters. This function
# uses breakpoint detection function writen by Dainius Masulinus and preprocessing + writen csv file by Sven-Arne

#' @param breaks number of breaks that are set of the algorithm e.g. BIC, LWZ, and RSS
#' @param magthreshold  threshold above which the breaks are kept. Breaks lower than the threshold are discarded
#' @param formula Formula passed to sctest() and bfastpp()
#' @param magcomponent components (possibly multiple for which to calculate the magnitudes (trend/harmonsin1/harmoncos3/...)
#' @param coefcomponent component (one!) for coefficient difference calculation
#' @param validating_data_path path to data need to be detected. This data is a time series and has a gpkg format
#' @param output_path path and name of output needed to be writen into CSV
#' @param cl number of cores used to run the data on your computers. It is important to setup paralell computation
#' @param plot Whether to call plot.bfast0n() on the output
#'
#' @return a data frame that contains output of bfast lite joined with their corresponding sample_ids
#' or NA if not enough observations/error in running the function.
cal_BFAST <- function( breaks, magthreshold, formula, magcomponent, coefcomponent, validating_data_path, output_path,plot, cl=mycores) {
  params <- list(
    preprocessing = list(
      interpolate = FALSE,
      seasonality = ""),
    bfastLiteInputs = list(
      InputTS="",
      scrange=NULL, 
      scsig=0.05, 
      breaks= breaks,
      sctype="OLS-MOSUM", 
      maginterval=0.1, 
      magcomponent=magcomponent,
      magstat="RMSD", 
      magthreshold= magthreshold, 
      coefcomponent=coefcomponent,
      coefthresholds=c(0, 0), 
      plot=plot, 
      quiet=TRUE, 
      order=3,
      formula=formula, 
      TargetYears=NULL,
      seasonfreq=0.5, 
      breaknumthreshold=Inf, 
      altformula=NULL
    )
  )
  # load an preprocess the input datafile 
  l8ndvi <- prepL8NDVIdataframe(validating_data_path)
  
  # prepare time series ------------------------------------------------------
  
  # get frequency table of acquisitions per year
  tab <- getFreqTabByYear(validating_data_path)
  
  # set the seasonality of the input time series
  
  if(!is.numeric(params$preprocessing$seasonality)){
    # define the time series frequency by the minimum number of aquisitions in the time series
    params$preprocessing$seasonality <- min(tab$Frequency[2:(length(tab$Frequency)-1)])
  }
  
  # make list with timeseries (ts) objects 
  print("Prepare time series list")
  tslist <- pblapply(1:nrow(l8ndvi), 
                     FUN = getTrimmedts, 
                     dframe = l8ndvi, 
                     lookuptable = tab, 
                     tsfreq = params$preprocessing$seasonality, 
                     interpolate = params$preprocessing$interpolate, 
                     cl = mycores)
  # --------- RUN BFAST Lite -------------------------------------------
  runBFASTLite <- function(i,timeserieslist, parameters){
    parameters$InputTS <- timeserieslist[[i]]
    return(do.call(UseBFASTLite, parameters))
  }
  
  bfastres <- pblapply(1:length(tslist),
                       FUN = runBFASTLite,
                       timeserieslist = tslist,
                       parameters = params$bfastLiteInputs,
                       cl = mycores)
  # join sample ids and centroid coordinates to the BFAST Monitor output
  print("Join output to sample ID information")
  theoutput <- pblapply(1:length(bfastres), 
                        FUN = joinOutput, 
                        algo_list = bfastres, 
                        dfids = l8ndvi, 
                        cl = mycores)
  
  # concatenate all 1-row data.frames into one large data.frame
  theoutput <- data.table::rbindlist(theoutput, fill = TRUE)
  
  # add prediction ids to each prediction. 
  theoutput <- theoutput %>% 
    group_by(sample_id) %>% 
    mutate(pred_id = row_number()*10) %>% 
    as.data.frame()
  
  
  
  #' if there is a magnitude to report, 
  #' fill the Magnitude column with "Mangitude.before"
  if ("Magnitude.before" %in% colnames(theoutput)) {
    theoutput <- theoutput %>% 
      mutate(Magnitude = if_else(!is.na(Magnitude.RMSD),Magnitude.RMSD, 0))}
  if ("tile" %in% colnames(validating_data_path)) {return(theoutput)} else {
    # make an output file name 
    outname <- ouputpath
    
    
    # export
    data.table::fwrite(theoutput,
                       file = outname,
                       showProgress = TRUE)
    return(theoutput)
  }
  
}
# This function used to write the output after running the dataset on BFAST lite. This fuction only
writeoutput <- function(outputbfast, ouputpath) {
  # join sample ids and centroid coordinates to the BFAST Monitor output
  print("Join output to sample ID information")
  theoutput <- pblapply(1:length(outputbfast), 
                        FUN = joinOutput, 
                        algo_list = outputbfast, 
                        dfids = l8ndvi, 
                        cl = mycores)
  
  # concatenate all 1-row data.frames into one large data.frame
  theoutput <- data.table::rbindlist(theoutput, fill = TRUE)
  
  # add prediction ids to each prediction. 
  theoutput <- theoutput %>% 
    group_by(sample_id) %>% 
    mutate(pred_id = row_number()*10) %>% 
    as.data.frame()
  
  
  
  #' if there is a magnitude to report, 
  #' fill the Magnitude column with "Mangitude.before"
  if ("Magnitude.before" %in% colnames(theoutput)) {
    theoutput <- theoutput %>% 
      mutate(Magnitude = if_else(!is.na(Magnitude.RMSD),Magnitude.RMSD, 0))}
  # make an output file name 
  outname <- ouputpath
  
  print(theoutput)
  # export
  data.table::fwrite(theoutput,
                     file = outname,
                     showProgress = TRUE)
}