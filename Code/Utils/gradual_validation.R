# ---- Scripting preparations ----------------------------------------------------------- 
# load packages 
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(data.table, tidyverse, pbapply)

# source Utils 
source(here("GitHub_code", "Utils.R"))



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