# ----  General info --------------------------
#' 
#' Name:             Duyen Bui
#' Study programs:   Master of Geo-information Science 
#' Student number:   1044767
#' 
#' Project:          Thesis: Compare the performance of multivariate and univariate time series models in detecting global land cover change
#' Supervisors:      Dainius Masiliunas & Nandin-Erdene Tsendbazar
#'
#' Script to calculate number TP, TN, FP, and FN for NDVI, and Satellite image dataset. The metrics such as F1 score, precision, and sensitivity are caculated.
#'
#'_____________________________________________________________________
# load packages

if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(here, devtools, sf, data.table, tidyverse, pbapply, parallel, data.table, hash)

# source external functions
setwd("/home/duyen/Master_Thesis/")
here::i_am("Code/Utils/gradual_validation.R")
debugSource("/home/duyen/Master_Thesis/Code/Utils/gradual_validation.R")
debugSource("/home/duyen/Master_Thesis/GitHub_code/BFASTutils.R")
# add progress bar option to show it in the terminal 
pbo <- pboptions(type="timer")
mycores <- detectCores()

# read all csv files that are the results of BFAST lite for calibration subdataset
NDVI_breakpoint_list <- list.files(path = "./cali_ouput/", pattern = "^_output_BFASTlite") 
cal_comb_ref <-read.csv("./Intermediate product/cal_compl_data.csv")
id_cal <- read.csv("./Intermediate product/cloud_free_product/DRMAT_final_sampleid_cal.csv")
cal_comb_ref_remain <- cal_comb_ref[cal_comb_ref$sample_id %in% id_cal$sample_id,]
# calculate confusion matrix for all calibration scenario. Put all the results in a list
list_breakpoint_output <- data.frame()
for (file in NDVI_breakpoint_list) {
  read_file <- read.csv(paste0("/home/duyen/Master_Thesis/cali_ouput/",file))
  validation <- validateAlgorithmTotal(ref_df = cal_comb_ref_remain, 
                                       algo_df = read_file, cl= mycores)
  stats <- myFPStats(validation, NewAccuracy = TRUE)
  list_breakpoint_output <- rbind(list_breakpoint_output,stats)
} 
# Write the list result into a data frame
# trim unnescessary part of the names of the files to create a new list of names for all the scenario of the data frame
# remove .csv from the file names
name_file <- sub(".csv", "", NDVI_breakpoint_list)
# select part of the name from _output_BFASTlite......_ (e.g _output_BFASTlite_BIC_SandT_025)
# and make them in a list
list_file <- list()
for (file in name_file) {
  lastchar <- nchar(file)
  file <- substring(file, 18, lastchar)
  list_file <- append(list_file, file)
}
# turn the list into a data frame and create a dataframe from all the scenario's of the calibration
file_name_df <- t(data.frame(list_file, row.names = NULL))
colnames(file_name_df) <- "Parameters"
breakpoint_stats_bfastlite <- cbind(file_name_df,list_breakpoint_output)
# write the dataframe result in a csv file
write.csv(breakpoint_stats_bfastlite, file = "./Output/_Output_calibration_bfastlite.csv", row.names = FALSE)

## RESULT: trend+harmon, BIC, 0.25

################ NDVI VALIDATION DATASET #######################################
# calulate the confusion matrix for validation sub-data set
id_val <- read.csv("./Intermediate product/cloud_free_product/DRMAT_final_sampleid_val.csv")
colnames(id_val) <- c("id", "sample_id")
val_comb_ref <- read.csv("./Intermediate product/val_compl_data.csv")
val_comb_ref_remain <- val_comb_ref[val_comb_ref$sample_id %in% id_val$sample_id,]
NDVI_val <- read.csv("./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_SandT_025_val.csv")

NDVI_val <- validateAlgorithmTotal(ref_df = val_comb_ref_remain, 
                                   algo_df = NDVI_val, cl= mycores)

NDVI_stats_val <- myFPStats(NDVI_val, NewAccuracy = TRUE)
# write the result into a csv
write.csv(NDVI_stats_val, file = "./Output/_output_BFASTlite_NDVI_stats_val.csv", row.names = FALSE)

####################### METHOD 1 ###############################################
# Select the most dominant floor year that contain more than or equal to 3 breakpoints 
# from all bands. Then choose the highest magnitude breakpoints.
####################### SR CALIBRATION DATASET #################################
# NOTE: this method did not get an actual calibration parameter due to cutting corner. SO it ran on different thresholds and the threshold value that equals 4 gave the highest F1 score
# Calculate change threshold based on chi square distribution
# threshold <- qchisq(0.05, 6)
# 
# test_thes <- qt(0.90, 6)
# threshold <- qnorm(.70, mean=0, sd=6)
########## Calibrate for BFAST lite aggregation
# path of the output of BFAST lite for the satellite calibration sub-data set
SR_breakpoint <- "./Intermediate product/cloud_free_product/__output_BFASTlite_breakpoint_SR_cal.gpkg"
# read all the name in the layer
SRNames <- st_layers(SR_breakpoint)$name
# read all the layers of the satellite image
SR_breakpoints <- lapply(SRNames, function(name) st_read(SR_breakpoint, layer=name))
# Read the .rds file
SR_cal_datetime <- readRDS("./Intermediate product/cloud_free_product/__output_BFASTlite_datetime_SR_cal.rds")
# Agrregate the results of all bands
probability <- c(0.5, 0.55, 0.6, 0.63, 0.65, 0.69, 0.7,0.8,0.9,0.95,0.99)
window <- c(6, 12)

stats_list <- list()
for (prob in probability) { 
  for (win in window) {
    result <- aggregate_all_bands(SR_breakpoints, SR_cal_datetime, prob, win)
    
    # Convert the results from dictionary to dataframe. Remove all sample ids that are not detected as change
    key_bp <- keys(result)
    values_bp <- values(result)
    max.length <- max(sapply(values_bp,length)) # since values of each sample id have different length, make all lists have the same length
    values_bp <- lapply(lapply(values_bp, unlist), "length<-", max.length)
    key_bp <- lapply(key_bp, function(x) rep(x, max.length))
    key_bp <- as.data.frame(t(as.data.frame(key_bp)))
    values_bp <- as.data.frame(t(as.data.frame(values_bp)))
    values_bp <- values_bp %>% pivot_longer(cols = everything(), names_to = "id", values_to = "Breakpoint")
    key_bp <- key_bp %>% pivot_longer(cols = everything(), names_to = "id", values_to = "sample_id")
    Cold_bp <- cbind(values_bp, key_bp)
    Cold_bp <- Cold_bp[, c("sample_id", "Breakpoint")]
    Cold_bp_ <- Cold_bp[!Cold_bp$Breakpoint %in% c("not included"), ]
    Cold_bp <- na.omit(Cold_bp_)
    SR_breakpoint_nobreak <- SR_breakpoints[[1]][SR_breakpoints[[1]]$Breakpoint %in% c(-9999, NA), c("sample_id", "Breakpoint")]
    Cold_bp_complete <- rbind(SR_breakpoint_nobreak,Cold_bp )
    Cold_bp_complete$sample_id <- as.integer(Cold_bp_complete$sample_id)
    Cold_bp_complete$Breakpoint <- as.numeric(Cold_bp_complete$Breakpoint)
    
    validation_SR_cal <- validateAlgorithmTotal(ref_df = cal_comb_ref_remain, 
                                                algo_df = Cold_bp_complete, cl= mycores)
    stats_SR_cal <- myFPStats(validation_SR_cal, NewAccuracy = TRUE)
    stats_list[[paste0("prob_", prob, "_win_", win)]] <- stats_SR_cal
  }
  
}

OutFile_agg_SR_cal <- paste("./cali_ouput/", "_output_BFASTlite_agg_cali_SR_cal.rds", sep="")
# Save the overarching list into a file
saveRDS(stats_list, file = OutFile_agg_SR_cal)


# plot example
Cold_bp_complete_bp <- Cold_bp_complete[!Cold_bp_complete$Breakpoint %in% c(-9999),]

SR_GPKG = "./Intermediate product/cloud_free_product/BFAST_lite_SR_cal.gpkg"
# Read each name layer of the satellite data
SRNames = st_layers(SR_GPKG)$name
SR = lapply(SRNames, function(name) st_read(SR_GPKG, layer=name))

# load an preprocess the input datafile 
l8ndvi <- lapply(SR[2:7], function(x) prepL8NDVIdataframe(x))  # load an preprocess the input datafile 


# prepare time series ------------------------------------------------------

# get frequency table of acquisitions per year
tab <- getFreqTabByYear("./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg")

# set the seasonality of the input time series
if(!is.numeric(params$preprocessing$seasonality)){
  # define the time series frequency by the minimum number of aquisitions in the time series
  params$preprocessing$seasonality <- min(tab$Frequency[2:(length(tab$Frequency)-1)])
}

# make list with timeseries (ts) objects 
print("Prepare time series list")
tslist <- lapply(l8ndvi, function(x) pblapply(1:nrow(x), 
                                              FUN = getTrimmedts, 
                                              dframe = x, 
                                              lookuptable = tab, 
                                              tsfreq = params$preprocessing$seasonality, 
                                              interpolate = params$preprocessing$interpolate, 
                                              cl = mycores))
index <- which(SR[[2]]$sample_id ==1743583)

# 
# install.packages("ggfortify")
# library(ggplot2)
# install.packages("strucchange")
# library(ggfortify)
# library(strucchange)

# Plot a graph with time_data1
plot(tslist[[1]][[index]],                           
     ylim = c(10000, 30000),
     type = "l",
     col = 1)


abline(v = Cold_bp_complete_bp[Cold_bp_complete_bp$sample_id == 1743583, ]$Breakpoint, col = c("darkgreen", "blue"), 
       lty = c(1, 2), lwd = c(1, 3)) 
# Add line graphs of other two dataset
lines(tslist[[3]][[6106]],                             
      type = "l",
      col = 3)
lines(tslist[[4]][[6106]],                             
      type = "l",
      col = 4)
lines(tslist[[5]][[6106]],                             
      type = "l",
      col = 5)
lines(tslist[[6]][[6106]],                             
      type = "l",
      col = 6)
lines(tslist[[2]][[6106]],                             
      type = "l",
      col = 2)
changeid <- cal_comb_ref[!cal_comb_ref$change_yr %in% c(NA), ]
cal_comb_ref_remain[cal_comb_ref_remain$sample_id == 1743583,]
cal_comb_ref_remain[cal_comb_ref_remain$sample_id == 1021184,]
# # aggregate all satellite imagery bands
# aggegrate_SR_output_BFASTlite <- aggregatefun(unique(SR_breakpoints[[1]]$sample_id), SR_breakpoints,4)
# aggegrate_SR_output_BFASTlite <- unique(aggegrate_SR_output_BFASTlite)
# write.csv(aggegrate_SR_output_BFASTlite, file = "./Intermediate product/cloud_free_product/_output_BFASTlite_aggegrate_SR_output_cal.csv", row.names = FALSE)
# aggerate_SR_ouput_stats <- read.csv("./Intermediate product/cloud_free_product/_output_BFASTlite_aggegrate_SR_output_cal.csv")
# # calulate the confusion matrix for calibration sub-data set
# aggerate_SR_ouput_stats <- validateAlgorithmTotal(ref_df = cal_comb_ref, 
#                                                   algo_df = aggerate_SR_ouput_stats, cl= mycores)
# aggerate_SR_ouput_stats_ <- myFPStats(aggerate_SR_ouput_stats, NewAccuracy = TRUE)
###################### SR VALIDATION dataset ###################################
# path of the output of BFAST lite for the satellite calibration sub-data set
SR_breakpoint <- "./Intermediate product/cloud_free_product/__output_BFASTlite_breakpoint_SR_val.gpkg"
# read all the name in the layer
SRNames = st_layers(SR_breakpoint)$name
# read all the layers of the satellite image
SR_breakpoints <- lapply(SRNames, function(name) st_read(SR_breakpoint, layer=name))
# aggregate all satellite imagery bands
nrow(SR_breakpoints[[1]])
aggegrate_SR_output_BFASTlite <- aggregatefun(unique(SR_breakpoints[[1]]$sample_id), SR_breakpoints,4)
aggegrate_SR_output_BFASTlite <- unique(aggegrate_SR_output_BFASTlite) # if there is any duplicated rows, then remove
write.csv(aggegrate_SR_output_BFASTlite, file = "./Intermediate product/cloud_free_product/_output_BFASTlite_aggegrate_SR_output_val.csv", row.names = FALSE)
aggerate_SR_ouput_stats <- read.csv("./Intermediate product/cloud_free_product/_output_BFASTlite_aggegrate_SR_output_val.csv")
# calulate the confusion matrix for calibration sub-data set
aggerate_SR_ouput_stats <- validateAlgorithmTotal(ref_df = val_comb_ref, 
                                                  algo_df = aggerate_SR_ouput_stats, cl= mycores)
aggerate_SR_ouput_stats_ <- myFPStats(aggerate_SR_ouput_stats, NewAccuracy = TRUE)
write.csv("./Output/_output_BFASTlite_SR_aggregate_MaxMag_M1_val.csv")


# ####################### METHOD 2 ###############################################
# # Sum value of all band for each date. The output is one layer that contains summary values from all bands. Use this layer to run on the BFAST Lite.
# # There is no need to calibrate in this method, but I run the calibration subdata set anyway
# #calibration data
# SR_comb_BFASTLite_output <- read.csv("./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_T_025_SR_comb_M2_cal.csv")
# SR_comb_BFASTLite_output_stats <- validateAlgorithmTotal(ref_df = cal_comb_ref, 
#                                                          algo_df = SR_comb_BFASTLite_output, cl= mycores)
# SR_comb_BFASTLite_output_stats_ <- myFPStats(SR_comb_BFASTLite_output_stats, NewAccuracy = TRUE)
# write.csv(SR_comb_BFASTLite_output_stats_, "./Intermediate product/cloud_free_product/_output_BFASTLite_BIC_T_025_SR_comb_M2_cal_stats.csv")
# #validation data
# SR_comb_BFASTLite_output_val <- read.csv("./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_T_025_SR_comb_M2_val.csv")
# SR_comb_BFASTLite_output_stats_val <- validateAlgorithmTotal(ref_df = val_comb_ref, 
#                                                              algo_df = SR_comb_BFASTLite_output_val, cl= mycores)
# SR_comb_BFASTLite_output_stats_val_ <- myFPStats(SR_comb_BFASTLite_output_stats_val, NewAccuracy = TRUE)
# # write the result into a csv
# write.csv(SR_comb_BFASTLite_output_stats_val_, file = "./Output/_output_BFASTlite_BIC_T_025_SR_comb_M2_val.csv", row.names = FALSE)
