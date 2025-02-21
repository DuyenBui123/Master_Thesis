# ---- Scripting preparations ----------------------------------------------------------- 
# load packages 
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(data.table, tidyverse, pbapply, hash)

# source Utils 
source("/home/duyen/Master_Thesis/GitHub_code/Utils.R")



# 4) Did BFAST predict a break in the given year? TRUE/FALSE.
# BreakTimes: all breaks detected by BFAST.
# TargetYear: time at which we want to test.
# Threshold: how much fuzziness is allowed to still consider a break detected.
# This function needs to be run for every year and every point.
# If there was no break for that year, and BFAST predicted one, should return FALSE.
IsBreakInTargetYear = function(BreakTimes, TargetYear, threshold=1)
{
  #if (is.na(TargetYear))
  #    TargetYear = 2016 # If there is no break, we look at whether we predicted a break in 2016
  # TODO: Needs to be updated; previously, lack of break meant that the break time would be set to NA,
  # but now it's no longer the case, it's change_at_300m = FALSE and reference_year=<year>
  return(any(BreakTimes > TargetYear - threshold & BreakTimes < TargetYear+threshold))
}

# 4b) Vectorised version (takes a list of MODDetectBreaks() output and a column of target years)
# Returns a column of whether BFSAT predicted the break at that time or not.
VectorisedIsBreakInTargetYear = function(BreakList, threshold=0.5, TY=TargetYears)
{
  i = 1:length(BreakList)
  return(sapply(i, function(i){return(IsBreakInTargetYear(BreakList[[i]], TY[i], threshold=threshold))}))
}

#' 4c) A new version that takes all years as input and returns TP/FP/TN/FN.
#' It is more accurate, as a point may no longer be both a TP and FP.
#' @param PredictionDates Vector with predicted dates of breaks (year fractions)
#' @param TruthDates      Vector with reference dates of breaks (year fractions)
#' @param threshold       How much to allow the predictions to deviate while still considered correct (±years)
#' @param period          Start and end date of the period of interest (without threshold)
#' @return Vector with four values: TP, FP, TN, FN.
BreakConfusionStats = function(PredictionDates, TruthDates, threshold=0.5, period=c(2016, 2019))
{
  # Remove all predictions that fall outside of the period
  if (length(PredictionDates) > 0)
    PredictionDates = PredictionDates[PredictionDates > min(period) - threshold & PredictionDates < max(period) + threshold]
  
  PointDistance = function(point, points) any(abs(point - points) <= threshold+1e-13)
  
  # Algorithm: each category requires one rule to calculate.
  # 1) False positives:
  TruthCloseToPred = if (length(TruthDates) <= 0) {
    # If there is no true break, all predictions are far from truth
    # and all predictions are false positives
    rep(FALSE, length(PredictionDates))
  } else if (length(PredictionDates) <= 0) {
    # If there are no predictions, there are no false positives
    logical(0)
  } else {
    # Is there a true break within threshold of each predicted break?
    sapply(PredictionDates, PointDistance, TruthDates)
  }
  # Any predictions not close to a real break is a false positive
  FP = sum(!TruthCloseToPred)
  
  # 2) False negatives:
  PredCloseToTruth = if (length(PredictionDates) <= 0) {
    # If there is no prediction, all true breaks are far from predictions
    # and all true breaks are false negatives
    rep(FALSE, length(TruthDates))
  } else if (length(TruthDates) <= 0) {
    # If there are no true breaks, there are no false negatives
    logical(0)
  } else {
    # Is there a predicted break within threshold of each true break?
    sapply(TruthDates, PointDistance, PredictionDates)
  }
  # Any true break not close to a prediction is a false negative
  FN = sum(!PredCloseToTruth)
  
  # 3) True positives:
  # Any break that is close to another from both sets
  TP = min(sum(TruthCloseToPred), sum(PredCloseToTruth))
  
  # 4) True negatives:
  # The rest, given the total number of items we should have.
  # e.g. 4 - (sum of above)
  MaxItems = floor(max(period) - min(period) )#+ 2*threshold)
  TN = MaxItems - sum(FP, FN, TP)
  # NOTE: TN is never reliable!
  
  # Check whether the output makes sense
  Result = c(TP=TP, TN=TN, FP=FP, FN=FN)
  if (any(Result < 0)) {
    #warning(paste("Negative validation values:", toString(Result)))
    
    while (min(Result) < 0) {
      # Just clamp to 0
      MinStat = which.min(Result)
      Result[MinStat] = Result[MinStat]+1
      # Alternatively distribute the negativity over the most positive values
      #MaxStat = which.max(Result)
      #Result[MaxStat] = Result[MaxStat]-1
    }
    
  }
  #if (sum(Result) > MaxItems) {
  #    stop(paste("More validation items than possible:", toString(Result)))
  #}
  
  return(Result)
}
# 5) Get statistics from the prediction result.
# This is generic enough to handle multiple formats.
FPStats = function(predictions, truth = NULL, round=3)
{
  # If we get a data.frame, try to automagically determine what is predicted and what is true
  if (is.data.frame(predictions) && is.null(truth))
  {
    if ("bfast_guess" %in% names(predictions))
    {
      if ("change_at_100m" %in% names(predictions)) {
        truth = predictions$change_at_100m == "yes"
      } else if ("changed" %in% names(predictions)) {
        truth = predictions$changed
      }
      predictions = predictions$bfast_guess
    }
  }
  # New version with stats already in the DF
  if (all(c("TP", "TN", "FP", "TP") %in% names(predictions))) {
    TruePositiveCount = ifelse(all(is.na(predictions$TP)), NA, sum(predictions$TP, na.rm=TRUE))
    FalsePositiveCount = ifelse(all(is.na(predictions$FP)), NA, sum(predictions$FP, na.rm=TRUE))
    TrueNegativeCount = ifelse(all(is.na(predictions$TN)), NA, sum(predictions$TN, na.rm=TRUE))
    FalseNegativeCount = ifelse(all(is.na(predictions$FN)), NA, sum(predictions$FN, na.rm=TRUE))
  } else {
    # We predicted a break and it was a break
    TruePositiveCount = sum(predictions & truth, na.rm=TRUE)
    # We predicted a break but there wasn't one (overprediction)
    FalsePositiveCount = sum(predictions & !truth, na.rm=TRUE)
    # We predicted no break, and there were none
    TrueNegativeCount = sum(!predictions & !truth, na.rm=TRUE)
    # We predicted no break, but there was one (we missed it)
    FalseNegativeCount = sum(!predictions & truth, na.rm=TRUE)
  }
  # Percent of true positives out of all change
  Sensitivity = TruePositiveCount / (TruePositiveCount + FalseNegativeCount) # AKA Recall, Previously TruePositiveRate
  Specificity = TrueNegativeCount / (TrueNegativeCount + FalsePositiveCount)
  # Percent of false positive out of no change
  FalsePositiveRate = FalsePositiveCount / (TrueNegativeCount + FalsePositiveCount) # False positive rate or alpha or p-value or Type I Error
  PositiveProportion = TruePositiveCount / FalsePositiveCount
  PositiveLikelihood = Sensitivity / FalsePositiveRate # Likelihood Ratio for Positive Tests
  Precision = TruePositiveCount / (TruePositiveCount + FalsePositiveCount) # AKA positive predictive value
  Accuracy = (TruePositiveCount + TrueNegativeCount) /
    sum(TruePositiveCount, FalsePositiveCount, TrueNegativeCount, FalseNegativeCount)
  F1Score = 2 * (Precision * Sensitivity)/ (Precision + Sensitivity)
  Beta = abs(Precision - Sensitivity)
  return(data.frame(TruePositiveCount, FalsePositiveCount, TrueNegativeCount, FalseNegativeCount,
                    Sensitivity=round(Sensitivity, round), Specificity=round(Specificity, round),
                    Precision=round(Precision, round), F1Score=round(F1Score, round),
                    FalsePositiveRate=round(FalsePositiveRate, round),
                    PositiveProportion=round(PositiveProportion, round),
                    PositiveLikelihood=round(PositiveLikelihood, round), Accuracy=round(Accuracy, round),
                    Beta = round(Beta, round)))
}


#' Function to get validation statistics per algorithm. This function gets total 
#' scores per validation type (e.g.: 2 TP, 1 FN, 1 TN) in stead of per year. 
#' 
#' @return data.frame 
#' @param ref_df: data.frame with reference data for every sample_id and every year. 
#'               Assumed to have at least columns: sample_id, change_at_100m, and change_yr
#' @param algo_df: data.frame with outputs of an algorithm. 
#'                Assumed to have at least the following columns: 
#'                - sample_id
#'                - Breakpoint: time of a break in decimal year format  
#' @param threshold numeric. Offset in years in which a break is considered to be detected correctly or not.
#' @param NewAccuracy Whether or not to use the new accuracy calculation method (4c)
validateAlgorithmTotal <- function(ref_df, algo_df, threshold=1, quiet=FALSE, 
                                   NewAccuracy=TRUE, fpstats = TRUE, verbose=FALSE, 
                                   cl=NULL, ...)
{
  if (!NewAccuracy) {
    ref_df$bfast_guess <- NA
  } else {
    ref_df$TP <- NA
    ref_df$TN <- NA
    ref_df$FP <- NA
    ref_df$FN <- NA
  }
  
  if (!quiet) {
    pbi <- 0
    pb <- txtProgressBar(pbi, length(unique(ref_df$sample_id)), style = 3)
  }
  
  # Detect the breaks in a loop over unique points
  ProcessSampleTotal <- function(i)
  {
    # get a vector of break years
    
    BreakTimes <- algo_df[algo_df$sample_id == i,]$Breakpoint
    SampleChunk <- ref_df[ref_df$sample_id == i, ]
    
    if (all(is.na(BreakTimes))){
      # if the algorithm returned NA for whatever reason, 
      # (like errors, too many clouds, or too little data)
      # we cannot get statistics from that point, so it will be NA
      if(!NewAccuracy){
        SampleChunk <- SampleChunk[,c("sample_id", "change_yr", "change_at_100m", "bfast_guess")]
      }
      return(SampleChunk)
    }
    
    if (NewAccuracy) {
      TruthDates <- SampleChunk[SampleChunk$change_at_100m == TRUE,]$change_yr
      
      if (length(BreakTimes) == 1){      #
        if (BreakTimes == -9999){
          # the algorithm did not detect a break, 
          # so BreakTimes will be set to 0 for 
          # the BreakConfusionStats function
          BreakTimes <- numeric(0)
        }      
        
      }
      # Run the function to get a named vector with TP, TN, FP, FN
      Stats <- BreakConfusionStats(BreakTimes, TruthDates, threshold = threshold, 
                                   period = c(2016, 2019))
      # change the values of each column 
      SampleChunk <- SampleChunk %>% mutate(TP = Stats['TP'],TN = Stats['TN'],
                                            FP = Stats['FP'],FN = Stats['FN'])
      
    } else {
      for (year in as.data.frame(ref_df)[ref_df$sample_id == i,"change_yr"])
      {
        # the name "bfast_guess" was chosen for compatibility with the FPStats function 
        SampleChunk[SampleChunk$change_yr == year, "bfast_guess"] <- IsBreakInTargetYear(BreakTimes, year)
        SampleChunk <- SampleChunk[,c("sample_id", "change_yr", "change_at_100m", "bfast_guess")]
      }
    }
    return(SampleChunk)
  }
  print("Processing confusion statistics for every sample_id")
  ProcessedDF <- pblapply(unique(ref_df$sample_id), ProcessSampleTotal, cl=cl)
  
  ProcessedDF <- data.table::rbindlist(ProcessedDF) %>% as.data.frame()
  
}


myFPStats <- function(prediction, NewAccuracy = TRUE) {
  
  
  # For the new style of accuracy assessment, we don't need repeated years
  if (NewAccuracy){
    prediction <- prediction[!duplicated(prediction$sample_id),c("sample_id","TP","TN","FP","FN")]
  }
  
  # run validation statisics 
  
  prediction <- FPStats(prediction)
  
  
  # this changes the confusing name "bfast_guess" in the data.frame to "algo_guess"
  if("bfast_guess" %in% colnames(prediction)){
    prediction <- prediction %>% rename_all(recode, bfast_guess = "algo_guess")
  }
  return(prediction)
}

#' Function to calculate validation statistics for COLD
#' @param cal_id_i: A data frame containing sample ids used for detectionA list of lists. Each list is a band. Each band contains a list of observation id sample sites. Each sample id contains detected breakpoint information
#' data.frame with outputs of an algorithm. 
#'                Assumed to have at least the following columns: 
#'                - sample_id
#'                - Breakpoint: time of a break in decimal year format
#'                - Magnitude  
#' @param cal_comb_ref_cold a dataframe with reference data of each sample ids used for detection
#' @param cold_result a list of file pathes that contain results of COLD detection.
#' Assump to contain followings columns:
#'  + ID: id was notted during detection in COLD implementation to track sample ids later
#'  + Breakpoint: fractional year representing for breakpoints
#' @param prob_calibrated a list of string that represent probability change used for the calibration. 
#' This is optional, used for creating name of each calibration setting

#' @return a list of two result. "stats": validation statistics,"cm": a dataframe containing breakpoints and its sample ids

COLD_stats <- function(cal_id_i, cal_comb_ref_cold, cold_result, prob_calibrated) {
  stats_cold_list <- list()
  for (file.nr in 1:length(cold_result)) {
    # read file
    bp_cold <- read.delim(paste0(path = "./Data/Data_COLD/", cold_result[[file.nr]]),header = F)
    bp_cold <- bp_cold[2:nrow(bp_cold),]
    colnames(bp_cold) <- c("ID","Breakpoint")
    bp_cold$sample_id <- NA
    # Assign sample id of cal_id_i to the dataframe of bp_cold(COLD breakpoint result) based on shared ID
    unique_id_sampleid <- cal_id_i[cal_id_i$ID %in% unique(bp_cold$ID),]
    bp_cold_merged <- merge(bp_cold, unique_id_sampleid, by = "ID", all.x = TRUE)
    bp_cold_merged <- bp_cold_merged[, c(4, 2)]
    colnames(bp_cold_merged) <- c("sample_id", "Breakpoint")
    # prepare a breakpoint table for validation
    # select sample ids that are not included in the detected break results
    notincl_sample_id <- cal_id_i$sample_id[!cal_id_i$sample_id %in% bp_cold_merged$sample_id ]
    notincl_sample_id <- as.data.frame(notincl_sample_id)
    colnames(notincl_sample_id) <- "sample_id"
    notincl_sample_id$Breakpoint <- -999
    # merge both sample id that detected to have break and the ones that are not detected with breaks
    bp_all_cold <- rbind(bp_cold_merged, notincl_sample_id)
    bp_all_cold$sample_id <- as.integer(bp_all_cold$sample_id)
    bp_all_cold$Breakpoint <- as.numeric(bp_all_cold$Breakpoint)
    # Validate the result
    SR_cal_COLD <- validateAlgorithmTotal(ref_df = cal_comb_ref_cold, 
                                          algo_df = bp_all_cold, cl= mycores)
    # statistical value
    SR_cal_COLD_stats <- myFPStats(SR_cal_COLD, NewAccuracy = TRUE)
    stats_cold_list[[paste0("prob_", prob_calibrated[[file.nr]])]] <- SR_cal_COLD_stats
    if (length(prob_calibrated) > 1) {
      statsandconfmtrx <- list("stats" = stats_cold_list, "cm" = bp_all_cold)
      
      return(statsandconfmtrx)
    } else {statsandconfmtrx <- list("stats" = SR_cal_COLD_stats, "cm" = bp_all_cold)
    return(statsandconfmtrx)}
    
  }
  
}

# bfast_lite_local_maxima <- function(date_time, sampleid_SR_cal, SR, spec_th ) {
#   sample_id <-c ()
#   date_bp <-c()
# 
#   ## turn datime into datafram and merge mag
#   for (id in sampleid_SR_cal) {
#     
#     ts_index <- which(sampleid_SR_cal == id )
#     bp <- lapply(SR, function(x) {y <- x[x$sample_id == id,]$Breakpoint
#     return(y)
#     })
#     if(any(as.integer(bp[[1]]) == -9999) &  any(as.integer(bp[[2]]) == -9999) & any(as.integer(bp[[3]]) == -9999) &
#        any(as.integer(bp[[4]]) == -9999) & any(as.integer(bp[[5]]) == -9999) & any(as.integer(bp[[6]]) == -9999)) {
#       sample_id <- c(sample_id, id)
#       date_bp <- c(date_bp, -9999)
#     } else {
#       ts_df <- as.data.frame(date_time)
#       colnames(ts_df) <- "datetime"
#       
#       merge_mag <- merge(ts_df, SR[[1]][SR[[1]]$sample_id == id,], by = "datetime", all.x= TRUE)
#       mag <- merge_mag$Magnitude
#       mag[is.na(mag)] <- 0
#       mag_ts <- ts(mag)
#       
#       
#       merge_mag_2 <- merge(ts_df, SR[[2]][SR[[2]]$sample_id == id,], by = "datetime", all.x= TRUE)
#       mag2 <- merge_mag_2$Magnitude
#       mag2[is.na(mag2)] <- 0
#       mag_ts2 <- ts(mag2)
#       
#       merge_mag_3 <- merge(ts_df, SR[[3]][SR[[3]]$sample_id == id,], by = "datetime", all.x= TRUE)
#       mag3 <- merge_mag_3$Magnitude
#       mag3[is.na(mag3)] <- 0
#       mag_ts3 <- ts(mag3)
#       merge_mag_4 <- merge(ts_df, SR[[4]][SR[[4]]$sample_id == id,], by = "datetime", all.x= TRUE)
#       mag4 <- merge_mag_4$Magnitude
#       mag4[is.na(mag4)] <- 0
#       mag_ts4 <- ts(mag4)
#       
#       merge_mag_5 <- merge(ts_df, SR[[5]][SR[[5]]$sample_id == id,], by = "datetime", all.x= TRUE)
#       mag5 <- merge_mag_5$Magnitude
#       mag5[is.na(mag5)] <- 0
#       mag_ts5 <- ts(mag5)
#       
#       merge_mag_6 <- merge(ts_df, SR[[6]][SR[[6]]$sample_id == id,], by = "datetime", all.x= TRUE)
#       mag6 <- merge_mag_6$Magnitude
#       mag6[is.na(mag6)] <- 0
#       mag_ts6 <- ts(mag6)
#       
#       suma <- mag_ts+ mag_ts2+mag_ts3+mag_ts4+mag_ts5+mag_ts6
#       suma <- as.vector(suma)
#       length(suma)
#       ts_sum <- c()
# 
#       for (i in 1:length(suma)) {
#         if (i+6< length(suma)) {
#           sum <- sum(suma[i:(i+6-1)])
#           ts_sum <- c(ts_sum, sum)
#         } else{
#           j = length(suma) - i
#           sum <- sum(suma[i:(i+j-1)])
#           ts_sum <- c(ts_sum, sum)
#         }
#       }
#       maxima <- localMaxima( ts_sum #ts_sum # change
#       )
#       maxima_old <- maxima
#       dis <- 1
#       while (dis == 1) {
#         if (length(maxima) ==1 ) {
#           sample_id <- c(sample_id, id)
#           date_bp <-c(date_bp, decimal_date(as.Date(ts_df[maxima,])))
#           break
#           
#         }else{
#           differ <- maxima[2:length(maxima)] - maxima[1:length(maxima)-1]
#           
#           if (all(differ>=dist)) {
#             final_maxima <- c()
#             for (i in maxima) {
#               if (ts_sum[i] >spec_th) { #change
#                 final_maxima <- c(final_maxima, i)
#               }
#               
#             }
#             if(length(final_maxima) == 0){
#               sample_id <- c(sample_id, id)
#               date_bp <-c(date_bp, -9999)
#               break
#             } else {
#               for (i in final_maxima) {
#                 sample_id <- c(sample_id, id)
#                 date_bp <-c(date_bp, decimal_date(as.Date(ts_df[i,])))
#                 
#               }
#             }
#             break
#             
#           } else {
#             ind <- c()
#             for (i in seq_along(differ)) {
#               if (differ[i] <dist) {
#                 ind <- c(ind, i)
#                
#               }
#             }
#               maxima_new <- maxima[-c(ind, ind +1)]
#   
#             for( i in ind) {
#               max_value <- max(ts_sum[maxima[i]], ts_sum[maxima[i+1]]) #change
#               max_index <- which(ts_sum == max_value)[[1]]#change
#               maxima_new <- c(maxima_new,max_index )
#             }
#             maxima_new <- sort(unique(maxima_new))
#             if(length(maxima_new) == 1) {
#               final_maxima <- c()
#               for (i in maxima_new) {
#                 if (ts_sum[i] >spec_th) { #change
#                   final_maxima <- c(final_maxima, i)
#                 }
#               }
#               final_maxima <- sort(unique(final_maxima))
#               if(length(final_maxima) == 0){
#                 sample_id <- c(sample_id, id)
#                 date_bp <-c(date_bp, -9999)
#                 break
#               } else {
#                 for (i in final_maxima) {
#                   sample_id <- c(sample_id, id)
#                   date_bp <-c(date_bp, decimal_date(as.Date(ts_df[i,])))
#                   
#                 }
#                 break
#               }
#               
#             } else {
#             differ <- maxima_new[2:length(maxima_new)] - maxima_new[1:length(maxima_new)-1]
#             if (any(differ <dist)) {
#               maxima <- maxima_new
#               dis <- 1
#               
#             } else {
#               final_maxima <- c()
#               for (i in maxima_new) {
#                 if (ts_sum[i] >spec_th) { #change
#                   final_maxima <- c(final_maxima, i)
#                 }
#               }
#               final_maxima <- sort(unique(final_maxima))
#               if(length(final_maxima) == 0){
#                 sample_id <- c(sample_id, id)
#                 date_bp <-c(date_bp, -9999)
#                 break
#               } else {
#                 for (i in final_maxima) {
#                   sample_id <- c(sample_id, id)
#                   date_bp <-c(date_bp, decimal_date(as.Date(ts_df[i,])))
#                   
#                 }
#                 break
#               }
#               
#             }
#             }
#           }
#           
#       
#           
# 
#         }
#         
#       } # end while loop
#     }# end else
#   } # end for loop
#   result <- list("sample_id" = sample_id, "bp" = date_bp)
#   return(result)
# }
#' Function to calculate local maxima
#' @param x is a vector
#' @return a vector of local maxima
localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

#' Function to aggregate results of all BFAST Lite single band
#' @param SR: A list of lists. Each list is a band. Each band contains a list of observation id sample sites. Each sample id contains detected breakpoint information
#' data.frame with outputs of an algorithm. 
#'                Assumed to have at least the following columns: 
#'                - sample_id
#'                - Breakpoint: time of a break in decimal year format
#'                - Magnitude  
#' @param sampleid_SR_cal a vector of sample ids which users want to test 
#' @param spec_th a numeric. Determining a threshold to discard noise or small breaks
#' @param span a numeric. Use to determine the observation numbers used to smooth a time series in LOESS
#' @param win a numeric. A window slide with size win. Use to make a sum of win consecutive points
#' @param dist a numeric. A distance between two consecutive breakpoints can occur
#' @return a list of two result. The first result is sample ids, the second one is breakpoints of the corresponding sample ids

bfast_lite_local_maxima_loess <- function(date_time, SR, sampleid_SR_cal, spec_th, span, win, dist ) {
  sample_id <-c ()
  date_bp <-c()
  
  # Loop through each id
  for (id in sampleid_SR_cal) {
    
    # check whether any band of an id contains bp. If not return no break
    bp <- lapply(SR, function(x) {y <- x[x$sample_id == id,]$Breakpoint
    return(y)
    })
    if(any(as.integer(bp[[1]]) == -9999) &  any(as.integer(bp[[2]]) == -9999) & any(as.integer(bp[[3]]) == -9999) &
       any(as.integer(bp[[4]]) == -9999) & any(as.integer(bp[[5]]) == -9999) & any(as.integer(bp[[6]]) == -9999)) {
      sample_id <- c(sample_id, id)
      date_bp <- c(date_bp, -9999)
    } else {# any band contains breakpoints
      ts_df <- as.data.frame(date_time)
      colnames(ts_df) <- "datetime"
      # create time series for each band which contain 0 for no break and value of magnitude for breaks
      merge_mag <- merge(ts_df, SR[[1]][SR[[1]]$sample_id == id,], by = "datetime", all.x= TRUE)
      mag <- merge_mag$Magnitude
      mag[is.na(mag)] <- 0
      mag_ts <- ts(mag)
      
      
      merge_mag_2 <- merge(ts_df, SR[[2]][SR[[2]]$sample_id == id,], by = "datetime", all.x= TRUE)
      mag2 <- merge_mag_2$Magnitude
      mag2[is.na(mag2)] <- 0
      mag_ts2 <- ts(mag2)
      
      merge_mag_3 <- merge(ts_df, SR[[3]][SR[[3]]$sample_id == id,], by = "datetime", all.x= TRUE)
      mag3 <- merge_mag_3$Magnitude
      mag3[is.na(mag3)] <- 0
      mag_ts3 <- ts(mag3)
      merge_mag_4 <- merge(ts_df, SR[[4]][SR[[4]]$sample_id == id,], by = "datetime", all.x= TRUE)
      mag4 <- merge_mag_4$Magnitude
      mag4[is.na(mag4)] <- 0
      mag_ts4 <- ts(mag4)
      
      merge_mag_5 <- merge(ts_df, SR[[5]][SR[[5]]$sample_id == id,], by = "datetime", all.x= TRUE)
      mag5 <- merge_mag_5$Magnitude
      mag5[is.na(mag5)] <- 0
      mag_ts5 <- ts(mag5)
      
      merge_mag_6 <- merge(ts_df, SR[[6]][SR[[6]]$sample_id == id,], by = "datetime", all.x= TRUE)
      mag6 <- merge_mag_6$Magnitude
      mag6[is.na(mag6)] <- 0
      mag_ts6 <- ts(mag6)
      # Sum all 6 time series
      suma_6 <- mag_ts+ mag_ts2+mag_ts3+mag_ts4+mag_ts5+mag_ts6
      suma_df <- as.data.frame(suma_6)
      colnames(suma_df) <- c("breakpoint")
      comb_bp_date <- cbind(ts_df, suma_df)
      comb_bp_date$datetime <- decimal_date(as.Date(comb_bp_date$datetime))
      # Use loess to smooth the time series
      comb_bp_date.lo <- loess(comb_bp_date$breakpoint~datetime , comb_bp_date, span = span, degree = 2)
      
      suma_fit <- predict(comb_bp_date.lo, comb_bp_date$datetime, se = TRUE)
      suma <- as.vector(suma_fit$fit)
      # make a time series with a window sliding which sums magnitude within the window size. 
      #Make a buffer in case detected breakpoints from all bands are not the same place
      ts_sum <- c()
      #plot(comb_bp_date$datetime, suma)
      for (i in 1:length(suma_6)) {
        if (i+win< length(suma_6)) {
          sum <- sum(suma_6[i:(i+win-1)])
          ts_sum <- c(ts_sum, sum)
        } else{
          j = length(suma_6) - i
          sum <- sum(suma_6[i:(i+j-1)])
          ts_sum <- c(ts_sum, sum)
        }
      }
      # calculate the local maxima. Use the value resulted from loess
      maxima <- localMaxima( suma #ts_sum # change
      )
      # While loop to check whether the distances between breaks are satisfied the requirement or not
      # select breakpoint that larger than the threshold
      dis <- 1
      while (dis == 1) {
        if (length(maxima) ==1 ) {
          sample_id <- c(sample_id, id)
          date_bp <-c(date_bp, decimal_date(as.Date(ts_df[maxima,])))
          break
          
        } else if(length(maxima) == 0) {sample_id <- c(sample_id, id)
        date_bp <-c(date_bp, -9999)
        break} else{#calculate the distance between two consecutive breakpoints
          differ <- maxima[2:length(maxima)] - maxima[1:length(maxima)-1]
          
          if (all(differ>=dist)) {
            final_maxima <- c()
            for (i in maxima) {
              if (ts_sum[i] >spec_th) { #change
                final_maxima <- c(final_maxima, i)
              }
              
            }
            if(length(final_maxima) == 0){
              sample_id <- c(sample_id, id)
              date_bp <-c(date_bp, -9999)
              break
            } else {
              for (i in final_maxima) {
                sample_id <- c(sample_id, id)
                date_bp <-c(date_bp, decimal_date(as.Date(ts_df[i,])))
                
              }
            }
            break
            
          } else {
            ind <- c()
            for (i in seq_along(differ)) {
              if (differ[i] <dist) {
                ind <- c(ind, i)
                
              }
            }
            maxima_new <- maxima[-c(ind, ind +1)]
            
            for( i in ind) {
              max_value <- max(ts_sum[maxima[i]], ts_sum[maxima[i+1]]) #change
              max_index <- which(ts_sum == max_value)[[1]]#change
              maxima_new <- c(maxima_new,max_index )
            }
            maxima_new <- sort(unique(maxima_new))
            if(length(maxima_new) == 1) {
              final_maxima <- c()
              for (i in maxima_new) {
                if (ts_sum[i] >spec_th) { #change
                  final_maxima <- c(final_maxima, i)
                }
              }
              final_maxima <- sort(unique(final_maxima))
              if(length(final_maxima) == 0){
                sample_id <- c(sample_id, id)
                date_bp <-c(date_bp, -9999)
                break
              } else {
                for (i in final_maxima) {
                  sample_id <- c(sample_id, id)
                  date_bp <-c(date_bp, decimal_date(as.Date(ts_df[i,])))
                  
                }
                break
              }
              
            } else { # after the first filter, if there are still breakpoints that are too closed to each other, then continue to filter them out
              differ <- maxima_new[2:length(maxima_new)] - maxima_new[1:length(maxima_new)-1]
              if (any(differ <dist)) {
                maxima <- maxima_new
                dis <- 1
                
              } else {
                final_maxima <- c()
                for (i in maxima_new) {
                  if (ts_sum[i] >spec_th) { #change
                    final_maxima <- c(final_maxima, i)
                  }
                }
                final_maxima <- sort(unique(final_maxima))
                if(length(final_maxima) == 0){
                  sample_id <- c(sample_id, id)
                  date_bp <-c(date_bp, -9999)
                  break
                } else {
                  for (i in final_maxima) {
                    sample_id <- c(sample_id, id)
                    date_bp <-c(date_bp, decimal_date(as.Date(ts_df[i,])))
                    
                  }
                  break
                }
                
              }
            }
          }
          
          
          
          
        }
        
      } # end while loop
    }
  }
  result <- list("sample_id" = sample_id, "bp" = date_bp)
  return(result)
}


# bfast_lite_local_maxima_loess_bk <- function(date_time, SR, sampleid_SR_cal, spec_th, span ) {
#   sample_id <-c ()
#   date_bp <-c()
#   
#   ## turn datime into datafram and merge mag
#   for (id in sampleid_SR_cal) {
#     
#     ts_index <- which(sampleid_SR_cal == id )
#     bp <- lapply(SR, function(x) {y <- x[x$sample_id == id,]$Breakpoint
#     return(y)
#     })
#     if(any(as.integer(bp[[1]]) == -9999) &  any(as.integer(bp[[2]]) == -9999) & any(as.integer(bp[[3]]) == -9999) &
#        any(as.integer(bp[[4]]) == -9999) & any(as.integer(bp[[5]]) == -9999) & any(as.integer(bp[[6]]) == -9999)) {
#       sample_id <- c(sample_id, id)
#       date_bp <- c(date_bp, -9999)
#     } else {
#       ts_df <- as.data.frame(date_time)
#       colnames(ts_df) <- "datetime"
#       
#       merge_mag <- merge(ts_df, SR[[1]][SR[[1]]$sample_id == id,], by = "datetime", all.x= TRUE)
#       mag <- merge_mag$Magnitude
#       mag[is.na(mag)] <- 0
#       mag_ts <- ts(mag)
#       
#       
#       merge_mag_2 <- merge(ts_df, SR[[2]][SR[[2]]$sample_id == id,], by = "datetime", all.x= TRUE)
#       mag2 <- merge_mag_2$Magnitude
#       mag2[is.na(mag2)] <- 0
#       mag_ts2 <- ts(mag2)
#       
#       merge_mag_3 <- merge(ts_df, SR[[3]][SR[[3]]$sample_id == id,], by = "datetime", all.x= TRUE)
#       mag3 <- merge_mag_3$Magnitude
#       mag3[is.na(mag3)] <- 0
#       mag_ts3 <- ts(mag3)
#       merge_mag_4 <- merge(ts_df, SR[[4]][SR[[4]]$sample_id == id,], by = "datetime", all.x= TRUE)
#       mag4 <- merge_mag_4$Magnitude
#       mag4[is.na(mag4)] <- 0
#       mag_ts4 <- ts(mag4)
#       
#       merge_mag_5 <- merge(ts_df, SR[[5]][SR[[5]]$sample_id == id,], by = "datetime", all.x= TRUE)
#       mag5 <- merge_mag_5$Magnitude
#       mag5[is.na(mag5)] <- 0
#       mag_ts5 <- ts(mag5)
#       
#       merge_mag_6 <- merge(ts_df, SR[[6]][SR[[6]]$sample_id == id,], by = "datetime", all.x= TRUE)
#       mag6 <- merge_mag_6$Magnitude
#       mag6[is.na(mag6)] <- 0
#       mag_ts6 <- ts(mag6)
#       
#       suma_6 <- mag_ts+ mag_ts2+mag_ts3+mag_ts4+mag_ts5+mag_ts6
#       suma_df <- as.data.frame(suma_6)
#       colnames(suma_df) <- c("breakpoint")
#       comb_bp_date <- cbind(ts_df, suma_df)
#       comb_bp_date$datetime <- decimal_date(as.Date(comb_bp_date$datetime))
#       # define function that returns the SSE
#       # calcSSE <- function(x){
#       #   loessMod <- try(loess(comb_bp_date$breakpoint~datetime, data=comb_bp_date, span=x), silent=T)
#       #   res <- try(loessMod$residuals, silent=T)
#       #   if(class(res)!="try-error"){
#       #     if((sum(res, na.rm=T) > 0)){
#       #       sse <- sum(res^2)  
#       #     }
#       #   }else{
#       #     sse <- 99999
#       #   }
#       #   return(sse)
#       # }
#       # 
#       # # Run optim to find span that gives min SSE, starting at 0.5
#       # optim(par=c(0.5, 0.1, 0.2, 0.3), calcSSE, method="SANN")
#       comb_bp_date.lo <- loess(comb_bp_date$breakpoint~datetime , comb_bp_date, span = span, degree = 2)
#       
#       suma_fit <- predict(comb_bp_date.lo, comb_bp_date$datetime, se = TRUE)
#       suma <- as.vector(suma_fit$fit)
#       length(suma)
#       ts_sum <- c()
#       #plot(comb_bp_date$datetime, suma)
#       for (i in 1:length(suma)) {
#         if (i+6< length(suma)) {
#           sum <- sum(suma[i:(i+6-1)])
#           ts_sum <- c(ts_sum, sum)
#         } else{
#           j = length(suma) - i
#           sum <- sum(suma[i:(i+j-1)])
#           ts_sum <- c(ts_sum, sum)
#         }
#       }
#       maxima <- localMaxima( suma #ts_sum # change
#       )
#       maxima_old <- maxima
#       if (length(maxima) ==1 ) {
#         sample_id <- c(sample_id, id)
#         date_bp <-c(date_bp, decimal_date(as.Date(ts_df[maxima,])))
#         
#       }else{
#         differ <- maxima[2:length(maxima)] - maxima[1:length(maxima)-1]
#         
#         if (all(differ>=dist)) {
#           final_maxima <- c()
#           for (i in maxima) {
#             if (suma_6[i] >spec_th) { #change
#               final_maxima <- c(final_maxima, i)
#             }
#             
#           }
#           if(length(final_maxima) == 0){
#             sample_id <- c(sample_id, id)
#             date_bp <-c(date_bp, -9999)
#             break
#           } else {
#             for (i in final_maxima) {
#               sample_id <- c(sample_id, id)
#               date_bp <-c(date_bp, decimal_date(as.Date(ts_df[i,])))
#               
#             }
#           }
#           
#           
#         } else {
#           ind <- c()
#           for (i in seq_along(differ)) {
#             if (differ[i] <dist) {
#               ind <- c(ind, i)
#               
#             }
#           }
#           maxima_new <- maxima[-c(ind, ind +1)]
#           
#           for( i in ind) {
#             max_value <- max(suma_6[maxima[i]], suma_6[maxima[i+1]]) #change
#             if (max_value == 0) {
#               max_index <- maxima[i+1]
#             } else {
#               max_index <- which(suma_6 == max_value)#change
#             }
#             maxima_new <- c(maxima_new,max_index)
#             
#           }
#           maxima_new <- sort(maxima_new)
#           
#           final_maxima <- c()
#           for (i in maxima_new) {
#             if (suma_6[i] >spec_th) { #change
#               final_maxima <- c(final_maxima, i)
#             }
#           }
#           final_maxima <- sort(final_maxima)
#           if(length(final_maxima) == 0){
#             sample_id <- c(sample_id, id)
#             date_bp <-c(date_bp, -9999)
#             break
#           } else {
#             for (i in final_maxima) {
#               sample_id <- c(sample_id, id)
#               date_bp <-c(date_bp, decimal_date(as.Date(ts_df[i,])))
#               
#             }
#           }
#         }
#         
#       }
#     }
#   }
#   result <- list("sample_id" = sample_id, "bp" = date_bp)
#   return(result)
# }