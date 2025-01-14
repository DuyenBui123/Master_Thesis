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
# add progress bar option to show it in the terminal 
pbo <- pboptions(type="timer")
mycores <- detectCores()

# Read reference data for calibration dataset
cal_comb_ref <- read.csv("./Intermediate product/cal_compl_data.csv")
# Read interpolated data with stl plus
cal_ndvi_nonNan <- read.csv("./Data/Data_DRMAT/ndvi_cal.csv")
# Remove sample ids that do not belong to the dataset used to detect breakpoint in DRMAT
cal_comb_ref <- tibble::rowid_to_column(cal_comb_ref, "ID")
# read sample ids that were not interpolated
rm_id <- read.csv("./Data/Data_DRMAT/deleted_sampleid_set.csv")
colnames(rm_id) <- NULL
rm_id <- rm_id[ , 2]
rm_id <- as.list(rm_id)
cal_comb_ref_rm_part <- cal_comb_ref[!cal_comb_ref$sample_id %in% rm_id,]
# calibrate for SR_cal interpolated by stlplus
# calibrate DRMAT with BIC, AIC, and HQC, and ridgid = 0.001, 0.0005, 0.005, 0.01, 0.1
# calculate the statistics for confusion matrices
bp_BIC_0001 <- DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_BIC_0001.txt", "./Data/Data_DRMAT/bpcm_S_BIC_0001.txt")
bp_BIC_00005 <- DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_BIC_00005.txt", "./Data/Data_DRMAT/bpcm_S_BIC_00005.txt")
bp_BIC_0005 <- DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_BIC_0005.txt", "./Data/Data_DRMAT/bpcm_S_BIC_0005.txt")
bp_BIC_001 <- DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_BIC_001.txt", "./Data/Data_DRMAT/bpcm_S_BIC_001.txt")


bp_AIC_00005 <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_AIC_00005.txt", "./Data/Data_DRMAT/bpcm_S_AIC_00005.txt")
bp_AIC_0005 <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_AIC_0005.txt", "./Data/Data_DRMAT/bpcm_S_AIC_0005.txt")
bp_AIC_0001 <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_AIC_0001.txt", "./Data/Data_DRMAT/bpcm_S_AIC_0001.txt")
bp_AIC_001 <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_AIC_001.txt", "./Data/Data_DRMAT/bpcm_S_AIC_001.txt")
bp_AIC_01 <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_AIC_01.txt", "./Data/Data_DRMAT/bpcm_S_AIC_01.txt")

bp_HQC_00005 <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_HQC_00005.txt", "./Data/Data_DRMAT/bpcm_S_HQC_00005.txt")
bp_HQC_0005 <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_HQC_0005.txt", "./Data/Data_DRMAT/bpcm_S_HQC_0005.txt")
bp_HQC_0001 <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_HQC_0001.txt", "./Data/Data_DRMAT/bpcm_S_HQC_0001.txt")
bp_HQC_001 <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_HQC_001.txt", "./Data/Data_DRMAT/bpcm_S_HQC_001.txt")
bp_HQC_01 <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_HQC_01.txt", "./Data/Data_DRMAT/bpcm_S_HQC_01.txt")
# Make a list of results of the calibration
list_of_cal <- rbind( bp_BIC_00005 = bp_BIC_00005,bp_BIC_0001 = bp_BIC_0001, bp_BIC_0005 =  bp_BIC_0005,bp_BIC_001 = bp_BIC_001,
                      bp_AIC_00005 = bp_AIC_00005, bp_AIC_0005 = bp_AIC_0005, bp_AIC_0001 =bp_AIC_0001,
                      bp_AIC_001 = bp_AIC_001,  bp_HQC_00005 = bp_HQC_00005,
                      bp_HQC_0005 = bp_HQC_0005, bp_HQC_0001 =bp_HQC_0001, bp_HQC_001 =bp_HQC_001)

write.csv(list_of_cal, "./Intermediate product/cloud_free_product/_output_DRMAT_SR_cal.csv" )
# calibrate for ndvi_cal interpolated by stlplus
bp_AIC_00005_ndvi <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_AIC_00005_ndvi.txt", "./Data/Data_DRMAT/bpcm_S_AIC_00005_ndvi.txt")
bp_BIC_00005_ndvi <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_BIC_00005_ndvi.txt", "./Data/Data_DRMAT/bpcm_S_BIC_00005_ndvi.txt", cal_comb_ref_rm_part, cal_ndvi_nonNan)
bp_HQC_00005_ndvi <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_HQC_00005_ndvi.txt", "./Data/Data_DRMAT/bpcm_S_HQC_00005_ndvi.txt")
# Make a list of results of the calibration
list_of_cal_ndvi <- rbind(bp_HQC_00005_ndvi = bp_HQC_00005_ndvi, bp_AIC_00005_ndvi = bp_AIC_00005_ndvi, bp_BIC_00005_ndvi = bp_BIC_00005_ndvi )
write.csv(list_of_cal_ndvi, "./Intermediate product/cloud_free_product/_output_DRMAT_ndvi_cal.csv" )

# calibrate for ndvi cal interpolated by stlplus + linear interpolation
cal_ndvi_nonNan_all <- read.csv("./Data/Data_DRMAT/ndvi_all_cal.csv")
cal_id <- read.csv("./Data/Data_DRMAT/id_ndvi_all_cal.csv")
cal_id <- cal_id[,2]
cal_comb_ref_rm <- cal_comb_ref[cal_comb_ref$sample_id %in% cal_id,]

bp_BIC_00005_ndvi_all <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_BIC_00005_all_ndvi.txt", "./Data/Data_DRMAT/bpcm_S_BIC_00005_all_ndvi.txt", cal_comb_ref_rm, cal_ndvi_nonNan_all)
bp_AIC_00005_ndvi_all <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_AIC_00005_all_ndvi.txt", "./Data/Data_DRMAT/bpcm_S_AIC_00005_all_ndvi.txt",  cal_comb_ref_rm, cal_ndvi_nonNan_all)
bp_HQC_00005_ndvi_all <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_HQC_00005_all_ndvi.txt", "./Data/Data_DRMAT/bpcm_S_HQC_00005_all_ndvi.txt", cal_comb_ref_rm, cal_ndvi_nonNan_all)
list_of_cal_ndvi_all <- rbind(bp_BIC_00005_ndvi_all = bp_BIC_00005_ndvi_all, bp_AIC_00005_ndvi_all = bp_AIC_00005_ndvi_all, bp_HQC_00005_ndvi_all = bp_HQC_00005_ndvi_all )
# calibrate for SR cal interpolated by stlplus + linear interpolation
bp_BIC_00005_all <- DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_BIC_00005_all.txt", "./Data/Data_DRMAT/bpcm_S_BIC_00005_all.txt", cal_comb_ref_rm, cal_ndvi_nonNan_all)
bp_AIC_00005_all <- DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_AIC_00005_all.txt", "./Data/Data_DRMAT/bpcm_S_AIC_00005_all.txt", cal_comb_ref_rm, cal_ndvi_nonNan_all)
bp_HQC_00005_all <- DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_HQC_00005_all.txt", "./Data/Data_DRMAT/bpcm_S_HQC_00005_all.txt", cal_comb_ref_rm, cal_ndvi_nonNan_all)
list_of_cal_all <- rbind(bp_BIC_00005_all = bp_BIC_00005_all, bp_AIC_00005_all = bp_AIC_00005_all, bp_HQC_00005_all = bp_HQC_00005_all )


#################### Validation data ###########################################
# BIC, 0.0005 give the highest score for all cases
# Read reference data for validation dataset
val_comb_ref <-read.csv("./Intermediate product/val_compl_data.csv")
val_ndvi_nonNan_all <- read.csv("./Data/Data_DRMAT/ndvi_all_val.csv")
val_comb_ref <- tibble::rowid_to_column(val_comb_ref, "ID")
val_id <- read.csv("./Data/Data_DRMAT/id_ndvi_all_val.csv")
val_id <- val_id[,2]
val_comb_ref_rm <- val_comb_ref[val_comb_ref$sample_id %in% val_id,]
# ndvi interpolated by stlplus + linear interpolation
bp_BIC_00005_all_ndvi_val <- DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_BIC_00005_all_ndvi_val.txt", "./Data/Data_DRMAT/bpcm_S_BIC_00005_all_ndvi_val.txt", val_comb_ref_rm, val_ndvi_nonNan_all)
write.csv(bp_BIC_00005_all_val, "./Output/_output_DRMAT_BIC_00005_all_val.csv")
# SR interpolated by stlplus + linear interpolation
bp_BIC_00005_all_val <- DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_BIC_00005_all_val.txt", "./Data/Data_DRMAT/bpcm_S_BIC_00005_all_val.txt", val_comb_ref_rm, val_ndvi_nonNan_all)
write.csv(bp_BIC_00005_all_ndvi_val, "./Output/_output_DRMAT_BIC_00005_all_ndvi_val.csv")


# SR interpolated by stlplus
val_ndvi_nonNan <- read.csv("./Data/Data_DRMAT/ndvi_val.csv")
val_id_notall <- read.csv("./Data/Data_DRMAT/id_ndvi_val.csv")
val_id_notall <- val_id_notall[,2]
val_comb_ref_rm_notall <- val_comb_ref[val_comb_ref$sample_id %in% val_id_notall,]

bp_BIC_00005_val <- DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_BIC_00005_val.txt", "./Data/Data_DRMAT/bpcm_S_BIC_00005_val.txt", val_comb_ref_rm_notall, val_ndvi_nonNan, val_id_notall)
write.csv(bp_BIC_00005_val, "./Output/_output_DRMAT_BIC_00005_SR_val.csv")


# # ndvi interpolated by stlplus
bp_BIC_00005_ndvi_val <- DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_BIC_00005_ndvi_val.txt",
                                      "./Data/Data_DRMAT/bpcm_S_BIC_00005_ndvi_val.txt", val_comb_ref_rm_notall, val_ndvi_nonNan, val_id_notall)
write.csv(bp_BIC_00005_ndvi_val, "./Output/_output_DRMAT_BIC_00005_ndvi_val.csv")

