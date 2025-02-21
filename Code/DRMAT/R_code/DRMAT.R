# ----  General info --------------------------
#' 
#' Name:             Duyen Bui
#' Study programs:   Master of Geo-information Science 
#' Student number:   1044767
#' 
#' Project:          Thesis: Compare the performance of multivariate and univariate time series models in detecting global land cover change
#' Supervisors:      Dainius Masiliunas & Nandin-Erdene Tsendbazar
#'
#' Script to calibrate the parameters for DRMAT and calculate validation data
#' to calculate confusion matrices and statistic
#' NOTE: the formula has changed, and was added more parameters. A few code was not adapted with newly added parameters.
#'
#'_____________________________________________________________________
# Load package
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(here, devtools, sf, data.table, tidyverse, pbapply, parallel, data.table, caret, quantmod,ff,foreign,R.matlab)
# Load external sources
setwd("/home/duyen/Master_Thesis")
source("./Code/Utils/gradual_validation.R")
source("./GitHub_code/BFASTutils.R")
debugSource("./Code/Utils/DRMAT_utils.R")
debugSource(here("Code", "Utils", "plot_utils.R"))
# add progress bar option to show it in the terminal 
pbo <- pboptions(type="timer")
mycores <- detectCores()
#################### Calibration data ##########################################
# Read reference data for calibration dataset
cal_comb_ref <-read.csv("./Intermediate product/cal_compl_data.csv")
id_cal <- read.csv("./Intermediate product/cloud_free_product/DRMAT_final_sampleid_cal.csv")
id_cal <- as.data.frame(id_cal[,2])
colnames(id_cal) <- c("sample_id")
cal_comb_ref_remain <- cal_comb_ref[cal_comb_ref$sample_id %in% id_cal$sample_id,]

# calibrate for SR_cal interpolated by stlplus
# calibrate DRMAT with BIC, AIC, and HQC, and ridgid = 0.001, 0.0005, 0.005, 0.01, 0.1
# calculate the statistics for confusion matrices
bp_BIC_00005 <-  DRMAT_conmax("./Data/Data_DRMAT/DRMAT_T_BIC_00005_SR_cal.txt", "./Data/Data_DRMAT/DRMAT_S_BIC_00005_SR_cal.txt", cal_comb_ref_remain, id_cal)
bp_BIC_0001 <-  DRMAT_conmax("./Data/Data_DRMAT/DRMAT_T_BIC_0001_SR_cal.txt", "./Data/Data_DRMAT/DRMAT_S_BIC_0001_SR_cal.txt", cal_comb_ref_remain, id_cal)
bp_BIC_001 <-  DRMAT_conmax("./Data/Data_DRMAT/DRMAT_T_BIC_001_SR_cal.txt", "./Data/Data_DRMAT/DRMAT_S_BIC_001_SR_cal.txt", cal_comb_ref_remain, id_cal)
ndvi_BIC_00005 <-  DRMAT_conmax("./Data/Data_DRMAT/DRMAT_T_BIC_00005_NDVI_cal.txt", "./Data/Data_DRMAT/DRMAT_S_BIC_00005_NDVI_cal.txt", cal_comb_ref_remain, id_cal)
ndvi_BIC_0001 <-  DRMAT_conmax("./Data/Data_DRMAT/DRMAT_T_BIC_0001_NDVI_cal.txt", "./Data/Data_DRMAT/DRMAT_S_BIC_0001_NDVI_cal.txt", cal_comb_ref_remain, id_cal)
ndvi_BIC_001 <-  DRMAT_conmax("./Data/Data_DRMAT/DRMAT_T_BIC_001_NDVI_cal.txt", "./Data/Data_DRMAT/DRMAT_S_BIC_001_NDVI_cal.txt", cal_comb_ref_remain, id_cal)



bp_AIC_00005 <-  DRMAT_conmax("./Data/Data_DRMAT/DRMAT_T_AIC_00005_SR_cal.txt", "./Data/Data_DRMAT/DRMAT_S_AIC_00005_SR_cal.txt", cal_comb_ref_remain, id_cal)
bp_AIC_0001 <-  DRMAT_conmax("./Data/Data_DRMAT/DRMAT_T_AIC_0001_SR_cal.txt", "./Data/Data_DRMAT/DRMAT_S_AIC_0001_SR_cal.txt", cal_comb_ref_remain, id_cal)
bp_AIC_001 <-  DRMAT_conmax("./Data/Data_DRMAT/DRMAT_T_AIC_001_SR_cal.txt", "./Data/Data_DRMAT/DRMAT_S_AIC_001_SR_cal.txt", cal_comb_ref_remain, id_cal)
ndvi_AIC_00005 <-  DRMAT_conmax("./Data/Data_DRMAT/DRMAT_T_AIC_00005_NDVI_cal.txt", "./Data/Data_DRMAT/DRMAT_S_AIC_00005_NDVI_cal.txt", cal_comb_ref_remain, id_cal)
ndvi_AIC_0001 <-  DRMAT_conmax("./Data/Data_DRMAT/DRMAT_T_AIC_0001_NDVI_cal.txt", "./Data/Data_DRMAT/DRMAT_S_AIC_0001_NDVI_cal.txt", cal_comb_ref_remain, id_cal)
ndvi_AIC_001 <-  DRMAT_conmax("./Data/Data_DRMAT/DRMAT_T_AIC_001_NDVI_cal.txt", "./Data/Data_DRMAT/DRMAT_S_AIC_001_NDVI_cal.txt", cal_comb_ref_remain, id_cal)

bp_HQC_00005 <-  DRMAT_conmax("./Data/Data_DRMAT/DRMAT_T_HQC_00005_SR_cal.txt", "./Data/Data_DRMAT/DRMAT_S_HQC_00005_SR_cal.txt", cal_comb_ref_remain, id_cal)
bp_HQC_0001 <-  DRMAT_conmax("./Data/Data_DRMAT/DRMAT_T_HQC_0001_SR_cal.txt", "./Data/Data_DRMAT/DRMAT_S_HQC_0001_SR_cal.txt", cal_comb_ref_remain, id_cal)
bp_HQC_001 <-  DRMAT_conmax("./Data/Data_DRMAT/DRMAT_T_HQC_001_SR_cal.txt", "./Data/Data_DRMAT/DRMAT_S_HQC_001_SR_cal.txt", cal_comb_ref_remain, id_cal)
ndvi_HQC_00005 <-  DRMAT_conmax("./Data/Data_DRMAT/DRMAT_T_HQC_00005_NDVI_cal.txt", "./Data/Data_DRMAT/DRMAT_S_HQC_00005_NDVI_cal.txt", cal_comb_ref_remain, id_cal)
ndvi_HQC_0001 <-  DRMAT_conmax("./Data/Data_DRMAT/DRMAT_T_HQC_0001_NDVI_cal.txt", "./Data/Data_DRMAT/DRMAT_S_HQC_0001_NDVI_cal.txt", cal_comb_ref_remain, id_cal)
ndvi_HQC_001 <-  DRMAT_conmax("./Data/Data_DRMAT/DRMAT_T_HQC_001_NDVI_cal.txt", "./Data/Data_DRMAT/DRMAT_S_HQC_001_NDVI_cal.txt", cal_comb_ref_remain, id_cal)
# Make a list of results of the calibration
list_of_cal <- rbind( bp_BIC_00005 = bp_BIC_00005,bp_BIC_0001 = bp_BIC_0001,bp_BIC_001 = bp_BIC_001,
                      bp_AIC_00005 = bp_AIC_00005, bp_AIC_0001 =bp_AIC_0001,
                      bp_AIC_001 = bp_AIC_001,  bp_HQC_00005 = bp_HQC_00005, bp_HQC_0001 =bp_HQC_0001, bp_HQC_001 =bp_HQC_001)

write.csv(list_of_cal, "./cali_ouput/_output_DRMAT_cali_SR_cal.csv" )
# Result calibration of SR cal: BIC:0.001

# Make a list of results of the calibration for ndvi
# Make a list of results of the calibration
list_of_ndvi_cal <- rbind( ndvi_BIC_00005 = ndvi_BIC_00005,ndvi_BIC_0001 = ndvi_BIC_0001,ndvi_BIC_001 = ndvi_BIC_001,
                           ndvi_AIC_00005 = ndvi_AIC_00005, ndvi_AIC_0001 =ndvi_AIC_0001,
                           ndvi_AIC_001 = ndvi_AIC_001,  ndvi_HQC_00005 = ndvi_HQC_00005, ndvi_HQC_0001 =ndvi_HQC_0001, ndvi_HQC_001 =ndvi_HQC_001)

write.csv(list_of_ndvi_cal, "./cali_ouput/_output_DRMAT_cali_ndvi_cal.csv" )
# Result calibration of ndvi: BIC 0.001

# # calibrate for ndvi cal interpolated by stlplus + linear interpolation
# cal_ndvi_nonNan_all <- read.csv("./Data/Data_DRMAT/ndvi_all_cal.csv")
# cal_id <- read.csv("./Data/Data_DRMAT/id_ndvi_all_cal.csv")
# cal_id <- cal_id[,2]
# cal_comb_ref_rm <- cal_comb_ref[cal_comb_ref$sample_id %in% cal_id,]
# 
# bp_BIC_00005_ndvi_all <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_BIC_00005_all_ndvi.txt", "./Data/Data_DRMAT/bpcm_S_BIC_00005_all_ndvi.txt", cal_comb_ref_rm, cal_ndvi_nonNan_all)
# bp_AIC_00005_ndvi_all <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_AIC_00005_all_ndvi.txt", "./Data/Data_DRMAT/bpcm_S_AIC_00005_all_ndvi.txt",  cal_comb_ref_rm, cal_ndvi_nonNan_all)
# bp_HQC_00005_ndvi_all <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_HQC_00005_all_ndvi.txt", "./Data/Data_DRMAT/bpcm_S_HQC_00005_all_ndvi.txt", cal_comb_ref_rm, cal_ndvi_nonNan_all)
# list_of_cal_ndvi_all <- rbind(bp_BIC_00005_ndvi_all = bp_BIC_00005_ndvi_all, bp_AIC_00005_ndvi_all = bp_AIC_00005_ndvi_all, bp_HQC_00005_ndvi_all = bp_HQC_00005_ndvi_all )
# # calibrate for SR cal interpolated by stlplus + linear interpolation
# bp_BIC_00005_all <- DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_BIC_00005_all.txt", "./Data/Data_DRMAT/bpcm_S_BIC_00005_all.txt", cal_comb_ref_rm, cal_ndvi_nonNan_all)
# bp_AIC_00005_all <- DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_AIC_00005_all.txt", "./Data/Data_DRMAT/bpcm_S_AIC_00005_all.txt", cal_comb_ref_rm, cal_ndvi_nonNan_all)
# bp_HQC_00005_all <- DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_HQC_00005_all.txt", "./Data/Data_DRMAT/bpcm_S_HQC_00005_all.txt", cal_comb_ref_rm, cal_ndvi_nonNan_all)
# list_of_cal_all <- rbind(bp_BIC_00005_all = bp_BIC_00005_all, bp_AIC_00005_all = bp_AIC_00005_all, bp_HQC_00005_all = bp_HQC_00005_all )


#################### Validation data ###########################################
# BIC, 0.0005 give the highest score for all cases
# Read reference data for validation dataset
val_comb_ref <-read.csv("./Intermediate product/val_compl_data.csv")
id_val <- read.csv("./Intermediate product/cloud_free_product/DRMAT_final_sampleid_val.csv")
id_val <- as.data.frame(id_val[,2])
colnames(id_val) <- c("sample_id")
val_comb_ref_remain <- val_comb_ref[val_comb_ref$sample_id %in% id_val$sample_id,]
# # ndvi interpolated by stlplus + linear interpolation
# bp_BIC_00005_all_ndvi_val <- DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_BIC_00005_all_ndvi_val.txt", "./Data/Data_DRMAT/bpcm_S_BIC_00005_all_ndvi_val.txt", val_comb_ref_rm, val_ndvi_nonNan_all)
# write.csv(bp_BIC_00005_all_val, "./Output/_output_DRMAT_BIC_00005_all_val.csv")
# # SR interpolated by stlplus + linear interpolation
# bp_BIC_00005_all_val <- DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_BIC_00005_all_val.txt", "./Data/Data_DRMAT/bpcm_S_BIC_00005_all_val.txt", val_comb_ref_rm, val_ndvi_nonNan_all)
# write.csv(bp_BIC_00005_all_ndvi_val, "./Output/_output_DRMAT_BIC_00005_all_ndvi_val.csv")


# # SR interpolated by stlplus
# val_ndvi_nonNan <- read.csv("./Data/Data_DRMAT/ndvi_val.csv")
# val_id_notall <- read.csv("./Data/Data_DRMAT/id_ndvi_val.csv")
# val_id_notall <- val_id_notall[,2]
# val_comb_ref_rm_notall <- val_comb_ref[val_comb_ref$sample_id %in% val_id_notall,]

bp_BIC_0001_val <- DRMAT_conmax("./Data/Data_DRMAT/DRMAT_T_BIC_0001_SRval.txt", "./Data/Data_DRMAT/DRMAT_S_BIC_0001_SRval.txt", val_comb_ref_remain, id_val)
write.csv(bp_BIC_0001_val, "./Output/DRMAT_SR_globe.csv", row.names = FALSE)
val_comb_ref_remain[val_comb_ref_remain$sample_id == 1405383,]
duyen <- bp_BIC_0001_val$cm[bp_BIC_0001_val$cm$sample_id == 1405383,]
bp_BIC_0001_val$stats
# # ndvi interpolated by stlplus
bp_BIC_0001_ndvi_val <- DRMAT_conmax("./Data/Data_DRMAT/DRMAT_T_BIC_0001_NDVI_val.txt",
                                     "./Data/Data_DRMAT/DRMAT_S_BIC_0001_NDVI_val.txt", val_comb_ref_remain, id_val)
bp_BIC_0001_ndvi_val$stats
write.csv(bp_BIC_0001_ndvi_val, "./Output/DRMAT_NDVI_globe.csv", row.names = FALSE)
############################## PLOT TIME SERIES ################################
# for SR
DRMAT_bp <- bp_BIC_0001_val$cm
DRMAT_bp_uni <- unique(DRMAT_bp$sample_id)
val_comb_ref_drmat_change <- val_comb_ref_remain[val_comb_ref_remain$label == "change",]
bp_id_uni <- unique(val_comb_ref_drmat_change$sample_id)
plot_ts_bp_act_cold_SR(1042255,DRMAT_bp, val_comb_ref_drmat_change, 
                       "./Intermediate product/cloud_free_product/DRMAT_SR_val.csv",
                       "./Intermediate product/cloud_free_product/DRMAT_NDVI_val.csv",
                       "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_val.gpkg")

DRMAT_ndvi_bp <- bp_BIC_0001_ndvi_val$cm
DRMAT_ndvi_bp_uni <- unique(DRMAT_ndvi_bp$sample_id)
plot_ts_bp_act_cold(1042255,DRMAT_ndvi_bp, val_comb_ref_drmat_change, 
                    
                    "./Intermediate product/cloud_free_product/DRMAT_NDVI_val.csv",
                    "./Intermediate product/cloud_free_product/BFAST_lite_NDVI_val.gpkg")

