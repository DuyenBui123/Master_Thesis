
# This function used to calculate confusion matrix and statistics of DRMAT output.
#' Breakpoints of seasonality and trend will be merged
#' 
#' @param bpTfile  a directory path of the breakpoints for trend signal as the output of DRMAT. It is supposed to be a text file
#' Suppose a dataframe includes a sample id and breakpoint columns which include fractional years for the breakpoint momenent 
#' @param bpSfile a directory path of the breakpoints for seasonality signal as the output of DRMAT. It is supposed to be a text file
#' Suppose a dataframe includes a sample id and breakpoint columns which include fractional years for the breakpoint momenent 
#' @param comb_ref_rm a reference dataset that only contains sample ids that used for detection
#' @param obs_nonNan observation dataset that was interpolated
#' @param sample_idsobs_nonNan a total reference dataset with all sample id
#' @return a dataframe result with statistics of confusion matrices
DRMAT_conmax <- function(bpTfile, bpSfile, comb_ref_rm, obs_nonNan, sample_ids_of_obs_nonNan) { # read the text file of seasonality breakpoints
  bp_S <- read.delim(bpSfile,header = F)
  # # read the text file of trend breakpoints
  bp_T <- read.delim(bpTfile,header = F)
  # Rearrange the table. A sample that has multiple breakpoint columns will be
  # converted to multiple rows (same sample id) with one columns (different breakpoints)
  bp_T_new <- bp_T %>% pivot_longer(!V1, names_to = "col", values_to = "breakpoint")
  # remove row with Nan
  bp_T_new <- na.omit(bp_T_new)
  # select only related columns
  bp_T_new <- subset(bp_T_new, select = -c(col))
  # change column names to match with the output of BFAST Lite
  colnames(bp_T_new) <- c("sample_id", "Breakpoint")
  # remove rows with -999
  bp_T_new <- bp_T_new[!bp_T_new$Breakpoint == -999,]
  # same procedure done for Seasonality
  bp_S_new <- bp_S %>% pivot_longer(!V1, names_to = "col", values_to = "breakpoint")
  bp_S_new <- subset(na.omit(bp_S_new), select = -c(col))
  colnames(bp_S_new) <- c("sample_id", "Breakpoint")
  bp_S_new <- bp_S_new[!bp_S_new$Breakpoint == -999,]
  # Combine trend and seasonality breakpoints
  bp_IC_ridgidnr <- rbind(bp_T_new, bp_S_new)
  # Add sample id that not contain breakpoints to complete the lists
  val_id <- as.data.frame(unique(sample_ids_of_obs_nonNan$sample_id)) # Select all sample_id from validation dataset
  bp_id <- as.data.frame((unique(bp_IC_ridgidnr$sample_id))) # Select all sample_id from the detected breakpoint list
  notadd_sampleid_bp <- as.data.frame(val_id[!val_id %in% bp_id]) # sample ids that not in the list will be saved
  colnames(notadd_sampleid_bp) <- c("sample_id")
  notadd_sampleid_bp$Breakpoint <- -999 # add -999 to the newly added sample ids as there are not breakpoints
  # colnames(notadd_sampleid_bp) <- c("sample_id", "Breakpoint")
  # bind beakpoint ids with ids no breakpoints
  bp_IC_ridgidnr_all <- rbind(bp_IC_ridgidnr,notadd_sampleid_bp)
  # validate the output with validation data
  DRMAT_SR_cal <- validateAlgorithmTotal(ref_df = comb_ref_rm, 
                                         algo_df = bp_IC_ridgidnr_all, cl= mycores)
  # calculate the stats
  DRMAT_IC_ridgidnr_stats <- myFPStats(DRMAT_SR_cal, NewAccuracy = TRUE)
  # return the stats
  return(DRMAT_IC_ridgidnr_stats)
  
  
}
# This function used to interpolate the missing data in the observation dataset using linear interpolation
#' 
#' @param rm_id  a directory path of the data that is needed to be interpolated. Suppose to be a time series data
#' @param tsorigin a ts list before interpolated but after preprocessed to determined seasonality parameter and remove dates 
#' if the number of dates per year exceeds the number of seasonality. In this case it is equal 22
#' @param tslist_notstpl a ts dataframe that contains all sample id that was not interpolated by stlplus. Need to 
#' have the same formate as ts origin
#' @param ndvi a logic (TRUE or FALSE). If it is a TRUE, it means that the rm_id list is not yet complete. Therefore need to extract sample id from the total rm_id list that can be stilled 
#' interpolated. It is a FALSE, the rm_id input is the result of ndvi dataset

#' @return a list of results include: time series before being interpolated but after preprocessing (have fewer dates than the original input),
#' time series after being interpolated, a list of index of sample id that is needed to be removed

linear_inter <- function(rm_id, tsorigin, tslist_notstpl, ndvi) {
  if (ndvi) {# create an empty list that collects indices that contain fewer than 12 consective missing data
    list_consec_nan_12 <- c()
    for (id in rm_id) { # loop through indices in rm_id (index that got removed from ndvi dataset)
      x = is.na(tsorigin[[id]]) # check whether there is na or not
      
      if (all(x) != TRUE) { # if not all Na, then calculate the frequency of each consecutive scenario
        x_1 <- rle(x)$lengths[rle(x)$values]
        X1_frq <- as.data.frame(table(x_1)) # convert into a frequency table
        colnames(X1_frq) <- c("consec_nan", "count") # change column names
        if (!any(as.integer(as.character(X1_frq$consec_nan)) >=12)) { # if there is any consecutive scenario that has more than or 
          # equal to 12 nans in a row, then do not collect that index. Otherwives, collect the index)
          
          
          list_consec_nan_12 <- c(list_consec_nan_12, id)} 
        
      }
    }
    
    list_4_lin_inter <- list() # create an empty list to collect sample ids that are interpolated linearly
    for (row in 1:length(tsorigin)) {# make a forloop to run through the whole ts origin 
      if (row %in% list_consec_nan_12) { # if the index of that row is contained in list_consec_nan_12, then interpolate 
        lin_inter <- na.approx(tsorigin[[row]], rule = 2, na.rm = FALSE) # linear interpolation
        list_4_lin_inter <- append(list_4_lin_inter, list(lin_inter)) # append
      }
    }
    list_4_lin_inter_df <- as.data.frame(do.call(rbind, list_4_lin_inter)) # convert a list into dataframe
    colnames(list_4_lin_inter_df) <- df$date # replace the column names with all dates in of time series
    tslist_non_and_lin <-  rbind(tslist_notstpl, list_4_lin_inter_df) # combine data that interpolated with stplus and linear interpolation
    mylist <- list("remain_id" = list_consec_nan_12, "tslist_all" = tslist_non_and_lin) # Make a list of results need to be return
    return(mylist)} else { # if ndvi = FALSE, used the rm_id output of ndvi dataset
      list_4_lin_inter <- list()
      for (row in 1:length(tsorigin)) {
        if (row %in% rm_id) {
          lin_inter <- na.approx(tsorigin[[row]], rule = 2, na.rm = FALSE)
          list_4_lin_inter <- append(list_4_lin_inter, list(lin_inter))
        }
      }
      list_4_lin_inter_df <- as.data.frame(do.call(rbind, list_4_lin_inter))
      colnames(list_4_lin_inter_df) <- df$date
      tslist_non_and_lin <-  rbind(tslist_notstpl, list_4_lin_inter_df)
      return(tslist_non_and_lin)
      
    }
  
}

# This function used to interpolate the missing data in the observation dataset using stlplus
#' 
#' @param datafile  a directory path of the data that is needed to be interpolated. Suppose to be a time series data

#' @return a list of results include: time series before being interpolated but after preprocessing (have fewer dates than the original input),
#' time series after being interpolated, a list of index of sample id that is needed to be removed

trans_multi_SR <-  function(datafile) {
  # set parameter of a time series object
  params <- list(
    preprocessing = list(
      interpolate = FALSE,
      seasonality = "")
  )
  # load an preprocess the input datafile 
  l8ndvi <- prepL8NDVIdataframe(datafile)
  
  # prepare time series ------------------------------------------------------
  
  # get frequency table of acquisitions per year
  tab <- getFreqTabByYear(datafile)
  
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
  # save time series to a new variable
  tslist_ori <- tslist
  # create empty lists to colect index that needed to be removed that have too 
  # many consecutive missing data in a row (>8 Nans). And a list of data that has
  # different length due to there are only Nans in the sample id.
  remov_index <- list()
  diff_length <- list()
  # loop through row by row of a time series (each row is a unique sample id)
  for (sample_id_indx in 1:length(tslist)) {
    sample_id_ts = tslist[[sample_id_indx]]
    # Parameters
    n <- 186      # Total number of dates for a time series
    n.p <- 22     # Seasonal period (e.g., weekly data in a year)
    
    
    # Generate cycle indices for the seasonal period. A ts will be segmented in n groups with size = n.p
    cycleSubIndices <- rep(1:n.p, ceiling(n / n.p))[1:n]
    # Check whether the ts dates are different than 186 ( a standard of nr dates after turn datafile to a ts object)
    if (length(sample_id_ts) != 186) {
      diff_length <- append(diff_length, sample_id_indx ) # if length is different, dave the sample id index
      next
    }
    # Check for the whole ts whether the consecutive missing data number is exceed the threshold (8).
    # if yes, collect the index
    if (any(by(sample_id_ts, list(cycleSubIndices), function(x) all(is.na(x))))) {remov_index <- append(remov_index,sample_id_indx )
    print("break")
    next # skip that sample id and continue futher
    
    }
    # Interpolate the remain sample id using stlplus. Decomposed to trend, seasonality, and error
    NDVI_stlpl <- stlplus(sample_id_ts, n.p = 22,
                          l.window = 23, t.window = 35, s.window = "periodic", s.degree = 1)
    # construct new values for the data by only sum trend and seasonality signals
    reconstruct_NDVI_stlpl <- seasonal(NDVI_stlpl) + trend(NDVI_stlpl)
    # Replace only NANs with interpolated values
    tslist[[sample_id_indx]] <- ifelse(is.na(sample_id_ts), reconstruct_NDVI_stlpl, sample_id_ts)
    
  }
  # save the interpolated ts into a new variable
  tslist_filled <- tslist
  # combine all index that cannot be interpolated
  remov_comb <- append(remov_index, diff_length)
  # Make a list of results that want to be returned
  mylist <- list("origin" = tslist_ori, "filled" = tslist_filled, "rm_id" = remov_comb)
  
  return(mylist)
  
}