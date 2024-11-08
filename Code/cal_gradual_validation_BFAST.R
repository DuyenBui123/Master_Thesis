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
p_load(here, devtools, sf, data.table, tidyverse, pbapply, parallel, data.table)

# source external functions
here::i_am("Code\\Utils\\gradual_validation.R")
source(here("Code","Utils", "gradual_validation.R"))
source(here("GitHub_code", "BFASTutils.R"))
# add progress bar option to show it in the terminal 
pbo <- pboptions(type="timer")
mycores <- detectCores()
setwd("C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\")
# read all csv files that are the results of BFAST lite for calibration subdataset
NDVI_breakpoint_list <- list.files(path = "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\", pattern = "^_output_BFASTlite") 
cal_comb_ref <-read.csv("C:\\Master_Thesis\\Intermediate product\\cal_compl_data.csv")
# calculate confusion matrix for all calibration scenario. Put all the results in a list
list_breakpoint_output <- data.frame()
for (file in NDVI_breakpoint_list) {
  read_file <- read.csv(file)
  validation <- validateAlgorithmTotal(ref_df = cal_comb_ref, 
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
write.csv(breakpoint_stats_bfastlite, file = "C:\\Master_Thesis\\Intermediate product\\_stats_breakpoints_bfastlite.csv", row.names = FALSE)



################ NDVI VALIDATION DATASET #######################################
# calulate the confusion matrix for validation sub-data set
val_comb_ref <-read.csv("C:\\Master_Thesis\\Intermediate product\\val_compl_data.csv")
NDVI_val <- read.csv("C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_T_025_val.csv")
NDVI_val <- validateAlgorithmTotal(ref_df = val_comb_ref, 
                                   algo_df = NDVI_val, cl= mycores)
NDVI_stats_val <- myFPStats(NDVI_val, NewAccuracy = TRUE)
# write the result into a csv
write.csv(NDVI_stats_val, file = "C:\\Master_Thesis\\Output\\_output_BFASTlite_NDVI_stats_val.csv", row.names = FALSE)

####################### SR CALIBRATION DATASET ##################################

# path of the output of BFAST lite for the satellite calibration sub-data set
SR_breakpoint <- "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\__output_BFASTlite_breakpoint_SR_cal.gpkg"
# read all the name in the layer
SRNames = st_layers(SR_breakpoint)$name
# read all the layers of the satellite image
SR_breakpoints <- lapply(SRNames, function(name) st_read(SR_breakpoint, layer=name))
# aggregate all satellite imagery bands
aggegrate_SR_output_BFASTlite <- aggregatefun(unique(SR_breakpoints[[1]]$sample_id), SR_breakpoints,4)
aggegrate_SR_output_BFASTlite <- unique(aggegrate_SR_output_BFASTlite)
write.csv(aggegrate_SR_output_BFASTlite, file = "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_aggegrate_SR_output_cal.csv", row.names = FALSE)
aggerate_SR_ouput_stats <- read.csv("C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_aggegrate_SR_output_cal.csv")
# calulate the confusion matrix for calibration sub-data set
aggerate_SR_ouput_stats <- validateAlgorithmTotal(ref_df = cal_comb_ref, 
                                                  algo_df = aggerate_SR_ouput_stats, cl= mycores)
aggerate_SR_ouput_stats_ <- myFPStats(aggerate_SR_ouput_stats, NewAccuracy = TRUE)


