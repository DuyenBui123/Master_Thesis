# ----  General info --------------------------
#' 
#' Name:             Duyen Bui
#' Study programs:   Master of Geo-information Science 
#' Student number:   1044767
#' 
#' Project:          Thesis: Compare the performance of multivariate and univariate time series models in detecting global land cover change
#' Supervisors:      Dainius Masiliunas & Nandin-Erdene Tsendbazar
#'
#' Script to validate the result of COLD
#'
#'_____________________________________________________________________
# Load package
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(here, devtools, sf, data.table, tidyverse, pbapply, parallel, data.table, caret, quantmod,ff,foreign,R.matlab)
# source external functions
setwd("/home/duyen/Master_Thesis/")
debugSource(here("GitHub_code", "Utils.R"))
debugSource(here("GitHub_code", "BFASTutils.R"))
source(here("GitHub_code", "cglops-change-detection", "src", "bfast-cal", "04-validation.r"))
debugSource(here("Code", "Utils", "gradual_validation.R"))
# add progress bar option to show it in the terminal 

pbo <- pboptions(type="timer")
mycores <- detectCores()
#################### Calibration data ##########################################
# prepare reference data for SR cal
# Read reference data for calibration dataset
cal_comb_ref <- read.csv("./Intermediate product/cal_compl_data.csv")
# Read a file contain the remain id after stlplus interpolation 
cal_id_i <- read.csv("./Data/Data_DRMAT/id_ndvi_cal.csv", header = F, sep = ",")
# read a file that contains the id list with Nan after using stlplus interpolation
rm_ID <- read.csv("./Data/Data_COLD/rm_id_w_nan_af_stlinter_SR_cal.csv")
# remove IDs that still contain NA after using stlplus interpolation
cal_id_i <- cal_id_i[!cal_id_i$V1 %in% rm_ID,]
cal_id_i<- cal_id_i[,"V2"]
cal_id_i <- as.data.frame(cal_id_i)
colnames(cal_id_i)<- c("sample_id")
# generate IDs for the remained obs after remove ID with nan
cal_id_i <- tibble::rowid_to_column(cal_id_i, "ID")
# Filter out ids in ref data that are not contain in id file 
cal_comb_ref_cold <- cal_comb_ref[cal_comb_ref$sample_id %in% cal_id_i,]


# read the result from COLD
# the results used default parameter
bp_cold <- read.delim("/home/duyen/Master_Thesis/Data/Data_COLD/breakpoint_COLD_SR_cal_default.txt",header = F)
bp_cold <- bp_cold[2:nrow(bp_cold),]
colnames(bp_cold) <- c("ID","Breakpoint")
bp_cold$sample_id <- NA
# merge sample id of cal_id_i to the dataframe of bp_cold based on shared ID
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
bp_all_cold$sample_id <- as.numeric(bp_all_cold$sample_id)
bp_all_cold$Breakpoint <- as.numeric(bp_all_cold$Breakpoint)
# Validate the result
SR_cal_COLD <- validateAlgorithmTotal(ref_df = cal_comb_ref_cold, 
                                      algo_df = bp_all_cold, cl= mycores)
# statistical value
SR_cal_COLD_stats <- myFPStats(SR_cal_COLD, NewAccuracy = TRUE)




