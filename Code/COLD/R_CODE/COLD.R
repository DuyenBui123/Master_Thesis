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
debugSource(here("Code", "Utils", "plot_utils.R"))
# add progress bar option to show it in the terminal 

pbo <- pboptions(type="timer")
mycores <- detectCores()
#################### Calibration data ##########################################
# prepare reference data for SR cal
# Read reference data for calibration dataset
cal_comb_ref <- read.csv("./Intermediate product/cal_compl_data.csv")
# Read a file contain the remain id after stlplus interpolation 
# sample id used for COLD detection
cal_id_i <- read.csv("./Intermediate product/cloud_free_product/DRMAT_final_sampleid_cal.csv", header = T, sep = ",")
cal_id_i <- as.data.frame(cal_id_i)
colnames(cal_id_i)<- c("ID", "sample_id")
# Filter out ids in ref data that are not contain in id file 
cal_comb_ref_cold <- cal_comb_ref[cal_comb_ref$sample_id %in% cal_id_i$sample_id,]
cal_comb_ref[cal_comb_ref$sample_id==1114775,]


# read the result from COLD
# the results used default parameter
# read all csv files that are the results of cold for calibration subdataset
COLD_breakpoint_list <- list.files(path = "./Data/Data_COLD/", pattern = "^breakpoint_COLD_SR_cal_") 
COLD_breakpoint_ndvi_list <- list.files(path = "./Data/Data_COLD/", pattern = "^breakpoint_COLD_ndvi_cal_")
prob_cal <- lapply(COLD_breakpoint_list, function(x) {return(substring(x, nchar(x) - 6))})
prob_ndvi_cal <- lapply(COLD_breakpoint_ndvi_list, function(x) {return(substring(x, nchar(x) - 6))})

stats_ndvi_cold <-COLD_stats(cal_id_i, cal_comb_ref_cold, COLD_breakpoint_ndvi_list,prob_ndvi_cal)
OutFile_COLD_ndvi_cal <- paste("./cali_ouput/", "_output_COLD_cali_ndvi_cal.rds", sep="")
# Save the overarching list into a file
saveRDS(stats_ndvi_cold, file = OutFile_COLD_ndvi_cal)
# result for ndvi calibration: prob: 0.90



stats_SRcal_cold <-COLD_stats(cal_id_i, cal_comb_ref_cold, COLD_breakpoint_list,prob_cal)
OutFile_COLD_SR_cal <- paste("./cali_ouput/", "_output_COLD_cali_SR_cal.rds", sep="")
# Save the overarching list into a file
saveRDS(stats_SRcal_cold, file = OutFile_COLD_SR_cal)
readRDS(paste("./cali_ouput/", "_output_COLD_cali_SR_cal.rds", sep=""))
# result for sr calibration: prob = 0.99

################## VALIDATION DATA
# prepare reference data for SR cal
# Read reference data for calibration dataset
val_comb_ref <- read.csv("./Intermediate product/val_compl_data.csv")
# Read a file contain the remain id after stlplus interpolation 
# sample id used for COLD detection
val_id_i <- read.csv("./Intermediate product/cloud_free_product/DRMAT_final_sampleid_val.csv", header = T, sep = ",")
val_id_i <- as.data.frame(val_id_i)
colnames(val_id_i)<- c("ID", "sample_id")
# Filter out ids in ref data that are not contain in id file 
val_comb_ref_cold <- val_comb_ref[val_comb_ref$sample_id %in% val_id_i$sample_id,]



# read the result from COLD
# the results used default parameter
# read all csv files that are the results of cold for calibration subdataset
COLD_breakpoint_list_val <- list.files(path = "./Data/Data_COLD/", pattern = "^breakpoint_COLD_srval_") 
COLD_breakpoint_ndvi_list_val <- list.files(path = "./Data/Data_COLD/", pattern = "^breakpoint_COLD_ndvival_")
prob_val <- lapply(COLD_breakpoint_list_val, function(x) {return(substring(x, nchar(x) - 6))})
prob_ndvi_val <- lapply(COLD_breakpoint_ndvi_list_val, function(x) {return(substring(x, nchar(x) - 6))})


COLD_breakpoint_ndvitest <- list.files(path = "./Data/Data_COLD/", pattern = "^breakpoint_COLD_ndvi_090_v1")
stats_ndvi_cold_val <- COLD_stats(val_id_i, val_comb_ref_cold, COLD_breakpoint_ndvitest,prob_ndvi_val)
write.csv(stats_ndvi_cold_val, "./Output/COLD_NDVI_globe.csv", row.names = FALSE)
stats_ndvi_cold_val$stats

stats_SRval_cold <- COLD_stats(val_id_i, val_comb_ref_cold, COLD_breakpoint_list_val,prob_val)
stats_SRval_cold$stats

write.csv(stats_SRval_cold, "./Output/COLD_SR_globe.csv", row.names = FALSE)
#################### PLOT ###########################
# for NDVI
COLD_bp <- stats_ndvi_cold_val$cm

val_comb_ref_cold_change <- val_comb_ref_cold[val_comb_ref_cold$label == "change",]
bp_id_uni <- unique(val_comb_ref_cold_change$sample_id)
COLD_bp_ndvi_uni <- unique(COLD_bp$sample_id)
plot_ts_bp_act_cold(1393825,COLD_bp, val_comb_ref_cold_change, 
                    "./Intermediate product/cloud_free_product/DRMAT_NDVI_val.csv", 
                    "./Intermediate product/cloud_free_product/BFAST_lite_NDVI_val.gpkg")
duyen <- val_comb_ref_cold[val_comb_ref_cold$sample_id %in% c(1367590, 1042255,1406712, 1367658, 1405383)                  ,]
write.csv(duyen, "./ref_sample_data.csv")
# for SR
val_comb_ref_cold[val_comb_ref_cold$sample_id == 1114775,]
COLD_SRbp <- stats_SRval_cold$cm
COLD_bp_SR_uni <- unique(COLD_SRbp$sample_id)
plot_ts_bp_act_cold_SR(1042255
                       ,COLD_SRbp, val_comb_ref_cold_change,
                       "./Intermediate product/cloud_free_product/DRMAT_SR_val.gpkg",
                       "./Intermediate product/cloud_free_product/DRMAT_NDVI_val.csv",
                       "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_val.gpkg")
# ########## test
# params <- list(
#   preprocessing = list(
#     interpolate = FALSE,
#     seasonality = ""),
#   bfastLiteInputs = list(
#     InputTS="",
#     scrange=NULL, 
#     scsig=0.05, 
#     breaks="LWZ",
#     sctype="OLS-MOSUM", 
#     maginterval=0.1, 
#     magcomponent= c("trend"),
#     magstat="RMSD", 
#     magthreshold=-Inf, 
#     coefcomponent=c("trend"),
#     coefthresholds=c(0, 0), 
#     plot=TRUE, 
#     quiet=TRUE, 
#     order=3,
#     formula=response ~ harmon + trend , 
#     TargetYears=NULL,
#     seasonfreq=3, 
#     breaknumthreshold=Inf, 
#     altformula=NULL
#   )
# )
# SR_GPKG = "./Intermediate product/cloud_free_product/DRMAT_SR_val.gpkg"
# # Read each name layer of the satellite data
# SRNames = st_layers(SR_GPKG)$name
# SR = lapply(SRNames, function(name) st_read(SR_GPKG, layer=name))
# val_ <- read.csv("./Intermediate product/cloud_free_product/DRMAT_NDVI_val.csv")
# l8ndvi <- lapply(SR[2:7], function(x) prepL8NDVIdataframe(x))
# 
# # prepare time series ------------------------------------------------------
# 
# # get frequency table of acquisitions per year
# tab <- getFreqTabByYear("./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_val.gpkg")
# 
# # set the seasonality of the input time series
# if(!is.numeric(params$preprocessing$seasonality)){
#   # define the time series frequency by the minimum number of aquisitions in the time series
#   params$preprocessing$seasonality <- min(tab$Frequency[2:(length(tab$Frequency)-1)])
# }
# 
# # make list with timeseries (ts) objects
# print("Prepare time series list")
# tslist <- lapply(l8ndvi, function(x) pblapply(1:nrow(x),
#                                               FUN = getTrimmedts,
#                                               dframe = x,
#                                               lookuptable = tab,
#                                               tsfreq = params$preprocessing$seasonality,
#                                               interpolate = params$preprocessing$interpolate,
#                                               cl = mycores))
