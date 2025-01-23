# load packages
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(here, lubridate)#, bfast, strucchangeRcpp)

#source utils 
source("/home/duyen/Master_Thesis/GitHub_code/Utils.R")
# library(strucchangeRcpp)
# library(bfast)
source("/home/duyen/Master_Thesis/GitHub_code/cglops-change-detection/src/utils/enable_fast_bfast.r")
source("/home/duyen/Master_Thesis/GitHub_code/cglops-change-detection/src/bfast-cal/plotting.r")
# debugSource(here("GitHub_code", "bfast", "R", "bfastpp.R"))
# debugSource(here("GitHub_code", "strucchangeRcpp", "R", "breakpoints.R"))
# debugSource(here("GitHub_code", "strucchangeRcpp", "R", "Fstats.R"))
# debugSource(here("GitHub_code", "strucchangeRcpp", "R", "RcppExports.R"))
# debugSource(here("GitHub_code", "strucchangeRcpp", "R", "critvals-monitoring.R"))
# debugSource(here("GitHub_code", "strucchangeRcpp", "R", "critvals.R"))
# debugSource(here("GitHub_code", "strucchangeRcpp", "R", "efp.R"))
# debugSource(here("GitHub_code", "strucchangeRcpp", "R", "gefp.R"))
# debugSource(here("GitHub_code", "strucchangeRcpp", "R", "matrix.R"))
# debugSource(here("GitHub_code", "strucchangeRcpp", "R", "monitoring.R"))
# debugSource(here("GitHub_code", "strucchangeRcpp", "R", "pvalue.Fstats.R"))
# debugSource(here("GitHub_code", "strucchangeRcpp", "R", "recresid.R"))
# debugSource(here("GitHub_code", "strucchangeRcpp", "R", "simul-monitoring.R"))
# debugSource(here("GitHub_code", "strucchangeRcpp", "R", "zzz.R"))
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
    bp <- try(strucchangeRcpp:::breakpoints(formula, data=bpp, h=h))
    # extract prediction value for the whole ts
    bp_pred <- strucchangeRcpp:::breakpoints_pred.formula(formula, data=bpp, h=h)
    if ("try-error" %in% class(bp)) {
      print("An error has occurred:")
      print(bp)
      return(Result <- data.frame("Breakpoint" = "Error", "Magnitude" = "Error", "Coefficients" = "Error", "Rsq" = "Error"))
    }
    
    #' now, obtain the timing of the breakpoints, 
    #' using the information criterion of `breaks`, which is `"LWZ"` by default
    bpOptim <- strucchangeRcpp:::breakpoints(bp, breaks=breaks) # Get breakpoint time
    
    # # adding rmse for each segment right before the breakpoints
    n <- bpOptim$nobs # number of observation after prepared by bpp func
    # check whether there is any break
    if(any(is.na(bpOptim$breakpoints))) { # if no break, only add the start and end of the ts
      nbp <- 0
      bp_rmse <- c(0, n)
    } else { # if there are breaks, add start and end obs to create complete segments of a ts
      nbp <- length(bpOptim$breakpoints)
      bp_rmse <- c(0, bpOptim$breakpoints, n)
    }
    # calculate adjusted root mean square error - predefined minimun rmse
    var_y <- bp$y[2:NROW(bp$y)] - bp$y[1:NROW(bp$y)-1]
    adj_rmse <- median(abs(var_y), 1)
    num_yrs <- 365.25 # determine number of days per year including leap years
    rmse <- c()
    
    if (length(bp_rmse) >2) { 
      for (nrb  in 1: (length(bp_rmse)-2)) { 
        segment <- bp_rmse[nrb]+1 : bp_rmse[nrb+1]-1
        # if the segment length is larger than 22 then calculate a temporal rmse
        if (length(segment) > 24) { segment_w_bp <- bp_rmse[nrb]+1 : bp_rmse[nrb+1]
        RSS_priorseg <- bp$RSS(bp_rmse[nrb]+1, bp_rmse[nrb+1]-1)
        datetime_seg <- bpp$datetime[segment_w_bp]
        datetime_seg <- gsub('-', '', datetime_seg)
        datetime_seg <- as.Date(datetime_seg, format="%Y%m%d")
        jul_date_seg <- yday(datetime_seg)
        d_rt <- jul_date_seg[1:NROW(jul_date_seg)] - jul_date_seg[NROW(jul_date_seg)]
        d_yr <- abs(round(d_rt/num_yrs)*num_yrs-d_rt)
        d_yr_sort <- order(d_yr[1:NROW(d_yr)-1])
        d_yr_selec <-  d_yr_sort[1:24]
        pred_seg <- bp_pred[d_yr_selec]
        act_seg <- bp$y[d_yr_selec]
        rmse_tem <- sqrt(sum((pred_seg - act_seg)^2))/sqrt(24-8)
        rmse <- c(rmse_tem,rmse)
        } else { # if the segment length is smaller than 22 then calculate a normal rmse
          
          RSS_priorseg <- bp$RSS(bp_rmse[nrb]+1, bp_rmse[nrb+1]-1)
          RMSE_priorseg <- sqrt(RSS_priorseg/length(segment))
          rmse <- c(RMSE_priorseg, rmse)
        }
        
      }
      # save date time for breakpoints
      bp_datetime <- c()
      # calculate change vector magnitude for each break
      vec_change <- c()
      for (break_nr in 1:(length(bp_rmse)-2)) {
        # get the minimun RMSE
        final_rmse <- max(rmse[break_nr], adj_rmse)
        vec_change <- c(abs(bp_pred[bpOptim$breakpoints[break_nr]] - bp$y[[bpOptim$breakpoints[break_nr]]])/final_rmse, vec_change)
        bp_datetime <- c(bpp$datetime[bpOptim$breakpoints[break_nr]], bp_datetime)
      }
    }
    
    
    
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
    if (length(vec_change) > 0) {
      # Convert the output into dates of break
      Result <- data.frame("Breakpoint" = bpp[bpOptim$breakpoints, "time"], 
                           "Magnitude" = bpMag, # Now all output of bpMag will be included in a new column
                           # and appended to the "Magnitude"-label automatically, 
                           # in case you want to output only 
                           # selected magnitude results, use column indexing, like: bpMag[,3]
                           "Coefficients" = bpCoef, "Rsq" = rsquared, "vecmag" = vec_change, "datetime" = bp_datetime) 
    } else {
      # Convert the output into dates of break
      Result <- data.frame("Breakpoint" = bpp[bpOptim$breakpoints, "time"], 
                           "Magnitude" = bpMag, # Now all output of bpMag will be included in a new column
                           # and appended to the "Magnitude"-label automatically, 
                           # in case you want to output only 
                           # selected magnitude results, use column indexing, like: bpMag[,3]
                           "Coefficients" = bpCoef, "Rsq" = rsquared, "vecmag" = NA, "datetime" =  NA) 
    }
    
    if (!is.null(bpMag) && !is.null(bpCoef)){
      # Keep breakpoints that are above the threshold of magnitude
      MagFilter <- abs(bpMag[,magstat]) > magthreshold
      # Keep breakpoints where the coefficient difference is big enough
      CoefFilter <- bpCoef < min(coefthresholds) | bpCoef > max(coefthresholds)
      Result <- Result[MagFilter & CoefFilter]
      Result$Rsq <- rsquared
      Result$vecmag <- vec_change
      
      #' The current code using the indexing of Result,
      #' might snap if more than one break is found.In that case, 
      #' try to select a column in the Magnitude output, like here below: 
      
      #data.frame("Breakpoint" = bpp[bpOptim$breakpoints, "time"][MagFilter & CoefFilter],
      #           "Magnitude" = bpMag[MagFilter][3],
      #           "Coefficients" = bpCoef[CoefFilter])
    }
    if (length(Result) < 1){
      return(Result <- data.frame("Breakpoint" = NA, "Magnitude" = NA, "Coefficients" = NA, "Rsq" = NA, "vegmag" = NA))} # If we filtered out all results, return FALSE again
    if (plot){ # For ones that got kept, plot black on top
      abline(v=Result, col="black")}
    return(Result)
    
  } else {
    if (!quiet){
      print("too cloudy")}
    return(Result <- data.frame("Breakpoint" = NA, "Magnitude" = NA, "Coefficients" = NA, "Rsq" = NA, "vegmag" = NA))
  }
}

UseBFASTLite_datetime <-  function(InputTS, scrange=NULL, scsig=0.05, breaks="LWZ",
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
    return(bpp$datetime)
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
  runBFASTLite_datetime <- function(i, timeserieslist, parameters) {
    
    parameters$InputTS <- timeserieslist[[i]]
    return(do.call(UseBFASTLite_datetime, parameters))
  }
  bfastres_datetime <- lapply(1:length(tslist), runBFASTLite_datetime, 
                              timeserieslist = tslist, 
                              parameters = params$bfastLiteInputs)
  
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
  if ("tile" %in% colnames(validating_data_path)) {
    bpanddatetime <- list("bp" = theoutput, "datetime" = bfastres_datetime)
    return(bpanddatetime)} else {
      # make an output file name 
      outname <- output_path
      
      
      # export
      data.table::fwrite(theoutput,
                         file = outname,
                         showProgress = TRUE)
      bpanddatetime <- list("bp" = theoutput, "datetime" = bfastres_datetime)
      return(bpanddatetime)
    }
  
}
# This function used to write the output after running the dataset on BFAST lite. This fuction only
writeoutput <- function(outputbfast, ouputpath) {
  
  outname <- ouputpath
  
  print(theoutput)
  # export
  data.table::fwrite(theoutput,
                     file = outname,
                     showProgress = TRUE)
}


#' Function to merge breakpoints from all bands together 

#' 
#' @return data.frame 
#' @param SR_breakpoints: a list of data.frame. Each list corresponds to a band (1,2,3,4,5,6). Each df contains sample ids and their corresponding breakpoints.
#' This df contains all sample ids which also do not contain breakpoints or could not be used to detected due to too much clouds. 
#' @param SR_cal_datetime : a list of bands. Each band contains lists of time series without NA. Each time series is corresponding with a sample id.
#' The order of the time series list matches with theirs corresponding sample id orders in SR_breakpoints   
#' @return a dictionary which contains sample ids and their corresponding breakpoints

aggregate_all_bands <- function(SR_breakpoints,SR_cal_datetime ) {
  sampleid_SR_cal <- unique(SR_breakpoints[[1]]$sample_id) # get the unique sample id
  SR_cal_datetimedict <- hash() # save a sample id and its ts in a dictionary
  merged_ts <- list()
  
  for (i in 1:length(sampleid_SR_cal)) {
    # Check if all the lists from list_2 to list_7 are identical for the i-th element
    if (identical(SR_cal_datetime$list_2[[i]], SR_cal_datetime$list_3[[i]]) && 
        identical(SR_cal_datetime$list_3[[i]], SR_cal_datetime$list_4[[i]]) && 
        identical(SR_cal_datetime$list_4[[i]], SR_cal_datetime$list_5[[i]]) && 
        identical(SR_cal_datetime$list_5[[i]], SR_cal_datetime$list_6[[i]]) &&
        identical(SR_cal_datetime$list_6[[i]], SR_cal_datetime$list_7[[i]])) {
      
      # if no difference, append any ts to the general list
      merged_ts <- append(merged_ts, list(SR_cal_datetime$list_2[[i]]))
      
      
    } else { # if not, extract the difference of other ts bands compared to the base band 2
      
      # Use lapply to find elements that are not in list_2 for each list
      extra_date <- lapply(SR_cal_datetime[-1], function(x) setdiff(x[[i]], SR_cal_datetime$list_2[[i]]))
      extra_date <- unique(extra_date)
      # Combine list_2's i-th element with all extra dates
      new_ts <- sort(c(SR_cal_datetime$list_2[[i]], unlist(extra_date)))
      # Add the new time series to merged_ts
      merged_ts <- append(merged_ts, list(new_ts))
    }
  }
  # merged_ts contain Nan and null data. Remove them
  deleted_indices <- which(sapply(merged_ts, is.null))
  merged_ts <- Filter(Negate(is.null), merged_ts)
  # remove the same id in sample id
  sampleid_SR_cal <- sampleid_SR_cal[-deleted_indices]
  # make dictionary for ts and sample ids
  for (nr.sampid in 1: length(sampleid_SR_cal)) {
    SR_cal_datetimedict[sampleid_SR_cal[nr.sampid]] <- merged_ts[[nr.sampid]]
  }
  # filter only sample id containing breaks
  SR_breakpoints_filled <- lapply(SR_breakpoints, function(band) band[!band$Breakpoint %in% c(-9999, NA),])
  
  final_bp <- hash()
  # loop through each sample id to aggregate the results from all band
  for (sample_id_nr in sampleid_SR_cal ) { print(paste("next sample id", sample_id_nr))
    bp_sp1_bands <- lapply(SR_breakpoints_filled, function(band) {
      band <- band[band$sample_id == sample_id_nr, c("sample_id", "vecmag", "datetime", "Breakpoint")]
      band$datetime <- as.Date(band$datetime)
      return(band)
    })
    
    SR_cal_datetime1 <- as.Date(unlist(values(SR_cal_datetimedict, keys = sample_id_nr)))
    
    
    # track whether there is any row in the sample id at all
    record <- 0
    if (nrow(bp_sp1_bands[[1]]) >0 | nrow(bp_sp1_bands[[2]]) >0|nrow(bp_sp1_bands[[3]]) >0|nrow(bp_sp1_bands[[4]]) >0|nrow(bp_sp1_bands[[5]]) >0|nrow(bp_sp1_bands[[6]]) >0) {
      # determine which band has the longest rows or breaks
      max_length <- which.max(c(NROW(bp_sp1_bands[[1]]), NROW(bp_sp1_bands[[2]]), NROW(bp_sp1_bands[[3]]), NROW(bp_sp1_bands[[4]]),
                                NROW(bp_sp1_bands[[5]]), NROW(bp_sp1_bands[[6]])))
      # loop through all bands till the longest breaks of a band gets to 0
      while (NROW(bp_sp1_bands[[max_length]] >= 1)) {selected_datetime <- hash()
      # select date time of the first break from all bands
      for (band in seq_along(bp_sp1_bands)) {
        if (length(bp_sp1_bands[[band]]$datetime) > 0) {
          selected_datetime[band] <- bp_sp1_bands[[band]]$datetime[1]
        }
      }
      
      record <- record + 1
      # Convert the date strings to Date objects
      dates <- as.Date(unlist(values(selected_datetime)))
      
      # Find the earliest date of the first break
      earliest_date <- which.min(dates)
      earliest_date_date <- min(dates, na.rm = TRUE)
      i <- 0
      # find the end date (+ 6 observations from the first break date)
      for (date. in SR_cal_datetime1) { i = i +1
      if ((earliest_date_date == as.Date(date.)) == TRUE) {
        end_date <- SR_cal_datetime1[i+6]
      }
      }
      
      # select dates that fall within the range dates
      selected_dates <- dates[dates >= earliest_date_date & dates <= end_date ]
      selected_dates_dict <- hash(selected_dates)
      
      # Get the keys (band number) corresponding to the selected dates
      selected_bands <-  as.numeric(unlist(keys(selected_dates_dict)))
      # calculate the change vector magnitude
      vec_mag <- c()
      for (i in selected_bands) {
        vec_mag <- c(bp_sp1_bands[[i]]$vecmag[1], vec_mag)
      }
      vec_mag_matrix <- matrix(vec_mag,nrow = 1,ncol = length(vec_mag))
      sum_vecmag <- norm(vec_mag_matrix)^2
      selected_vecmag <- hash()
      # detect break
      if (sum_vecmag > threshold) {
        for (band in seq_along(bp_sp1_bands)) {
          if (length(bp_sp1_bands[[band]]$datetime) > 0) {
            selected_vecmag[band] <- bp_sp1_bands[[band]]$vecmag[1]
          }
        }
        # Find the name (key) with the maximum value
        # Find the key with the maximum value
        max_key <- keys(selected_vecmag)[which.max(unlist(values(selected_vecmag)))]
        
        # Print the key with the maximum value
        max_key <- as.numeric(max_key)
        
        
        # save breakpoints that has largest vector magnitude
        if (all(has.key(as.character(bp_sp1_bands[[max_key]]$sample_id[1]), final_bp))== TRUE) {
          final_bp[[as.character(bp_sp1_bands[[max_key]]$sample_id[1])]] <- c(final_bp[[as.character(bp_sp1_bands[[max_key]]$sample_id[1])]], bp_sp1_bands[[max_key]]$Breakpoint[1])
        } else {final_bp[bp_sp1_bands[[max_key]]$sample_id[1]] <- bp_sp1_bands[[max_key]]$Breakpoint[1]
        }
        
      } else {final_bp[sample_id_nr] <- -9999 
      }
      # drop rows that used for calculating vec_mag
      for (bandnr. in selected_bands) {
        bp_sp1_bands[[bandnr.]] <- bp_sp1_bands[[bandnr.]][-1,]
      }
      } # end while loop
    } else {if (record == 0) { final_bp[sample_id_nr] <- "not included"
    next } else {
      next
    }}
  } # end for loop
  return(final_bp)
}