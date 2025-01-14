# ----  General info --------------------------
#' 
#' Name:             Duyen Bui
#' Study programs:   Master of Geo-information Science 
#' Student number:   1044767
#' 
#' Project:          Thesis: Compare the performance of multivariate and univariate time series models in detecting global land cover change
#' Supervisors:      Dainius Masiliunas & Nandin-Erdene Tsendbazar
#'
#' Script to interpolate missing data for ndvi and SR using stplus and linear interpolation
#' to convert data to the write format for DRMAT and COLD
#' to plot numbers of consecutive nan scenario
#'
#'_____________________________________________________________________
# load packages
# Load library 
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(abind, stlplus, TSstudio, zoo, here, devtools, sf, data.table, tidyverse, pbapply, parallel, data.table, ggplot2)
# Load external sources
source("/home/duyen/Master_Thesis/GitHub_code/BFASTutils.R")
source("/home/duyen/Master_Thesis/GitHub_code/getTrimmedts.R")
debugSource("/home/duyen/Master_Thesis/Code/Utils/Utils.R")
source("/home/duyen/Master_Thesis/Code/Utils/DRMAT_utils.R")
# Set working directory
setwd("/home/duyen/Master_Thesis/")
# create a new folfer
# setting up the main directory
main_dir <- "/home/duyen/Master_Thesis/Data/"
# setting up the sub directory
sub_dir <- "Data_DRMAT"
# check if sub directory exists 
if (file.exists(sub_dir)){
  
  # specifying the working directory
  setwd(file.path(main_dir, sub_dir))
} else {
  
  # create a new sub directory inside
  # the main path
  dir.create(file.path(main_dir, sub_dir))
  
  # specifying the working directory
  setwd(file.path(main_dir, sub_dir))
}
# Determine core for calculation
pbo <- pboptions(type="timer")
mycores <- detectCores()
###################### USING STLPLUS TO INTERPOLATE MISSING DATA FOR NDVI_CAL #################################

# Read ndvi files for calibration dataset with cloudfree data
l8ndvi <- st_read("./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg")
# Interpolate missing data using stlplus
ndvid_filled_cal <- trans_multi_SR("./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg")
# extract the results of the interpolated data: 1, a list of id that will not be interpolated and 2, ts with interpolated values
remov_comb_ndvi <- ndvid_filled_cal$rm_id
tslist_filled <- ndvid_filled_cal$filled # This ts list still include sample ids that have different length of columns and all Nans
# Plot time series
# test <- ts(tslist_spid2, start  = c(2013, 4), end = c(2021, 12), frequency = 22)
# ts.plot(test,tslist[[2]],
#         col = 2:1)

# These steps help to track which sample id index is interpolated and which ones need to be removed
# assign sample id again to the tslist that got interpolated
# Track sample ids which ones are remained
# turn ts lists into dataframe (result of interpolation)
tslist_filled_df <- do.call(rbind.data.frame, tslist_filled)
# bind sample id column to that ts df
tslist_filled_df <- cbind(tslist_filled_df, sample_id = l8ndvi$sample_id)
# create an id column for that ts df
tslist_filled_df$id <- seq.int(nrow(tslist_filled_df))

# remove the sample ids that could not be interpolated from the ts dataframe based in the id columns and rm_id output of the interpolation
# create a ts df without all uninterpolated sample ids
tslist_filled_df_rm <- tslist_filled_df[(!tslist_filled_df$id %in% remov_comb_ndvi),] # this df only used for extract interpolated sample_id


# Append rows that were interpolated in a new list and later use this list for create final observation data

for(row in 1:length(tslist_filled)) {if (length(tslist_filled[[row]]) == 186) {
  df <- as.data.frame(tslist_filled[[1]]) # extract colnames of the ts list
  df$date <- rownames(df)  # Add row names as a new column
  rownames(df) <- NULL     # Optional: Remove row names for cleaner data frame
  ts_notcomb <- df$date 
  break
}}

# assign filled ts list (all sample ids) to a new  variable
tslist_ndvi_186 <- tslist_filled
# create an empty list to collect interpolated sample ids
tslist_new <- list()
for (row in 1:length(tslist_ndvi_186)) { # only row index that not in the remove id list is selected
  # NOTE: this step seems to similar to the step done for the tslist_filled_df_rm. check
  if(!row %in% remov_comb_ndvi) { tslist_new <- append(tslist_new,list(tslist_ndvi_186[[row]]))}
}
# convert the final list into dataframe
tslist_new_df <- do.call(rbind.data.frame, tslist_new) # this is the final observation ts data 
# replace colnames with the dates of ts
colnames(tslist_new_df) <- df$date


# prepare sample id, time, ymd for DRMAT
# Use tslist_filled_df which contain remaining sample ide for creating a sample id dataframe
# created a dataframe of sample_ids that got interpolated with stlplus
sample_id_ndvi_cal <- get_sampleid(tslist_filled_df_rm, "Data_DRMAT/id_ndvi_cal.csv")
# create ymd of the time series dates (186)
ymd_ndvi_cal <- getymd(tslist_new_df, "Data_DRMAT/ymd_ndvi_cal.csv" )
# create a dataframe with a fractional year for each date in the time series dates
t_ndvi_cal <- getfracyear(tslist_new_df, "Data_DRMAT/t_ndvi_cal.csv")
# write the final observation data with interpolated values done by stlplus
outputname_ndvi_cal <- paste0(main_dir, "Data_DRMAT/ndvi_cal.csv")
write.csv(tslist_new_df, outputname_ndvi_cal)


######################## USING LINEAR INTERPOLATION IN NDVI_CAL DATA ###########
# interpolate data that was not interpolated by stlplus
tslist_non_and_lin_inter <- linear_inter(remov_comb_ndvi, ndvid_filled_cal$origin, tslist_new_df, ndvi = TRUE)


# adding ids to filled ts df (all)
if (nrow(tslist_filled_df) == 19251) {
  tslist_filled_df$id <- seq.int(nrow(tslist_filled_df))} else {
    print("Total number of the observation is not equal 19251")
  }
# get sample ids that were selected to be linearly interpolated
sample_id_lininter <- as.data.frame((tslist_filled_df[tslist_filled_df$id %in% tslist_non_and_lin_inter$remain_id, ]$sample_id))
colnames(sample_id_lininter) <- c("sample_id")
# sample ids of stlplus
colnames(sample_id_ndvi_cal) <- c("sample_id")
# append sample id of stlplus and linear interpolation
sample_id_ndvi_all_cal <- as.data.frame(c(as.integer(sample_id_ndvi_cal$sample_id), sample_id_lininter$sample_id))
colnames(sample_id_ndvi_all_cal) <- c("sample_id")
sample_id_ndvi_all_cal <- as.integer(sample_id_ndvi_all_cal$sample_id)
# write new sample id for stlplus+linear (all-except for all nans?) and new ndvi cal (all interpolated-except for all nans?)
write.csv(sample_id_ndvi_all_cal, paste0(main_dir, "Data_DRMAT/id_ndvi_all_cal.csv"))
write.csv(tslist_non_and_lin_inter$tslist_all, paste0(main_dir, "Data_DRMAT/ndvi_all_cal.csv"))

############################## USING STLPLUS + LINEAR TO INTERPOLATE MISSING DATA FOR SR CAL##########################
 # Read the CALIBRATION set of satellite data 
 SR_GPKG = "./Intermediate product/cloud_free_product/__cloudfree_L8TS__SR_cal.gpkg"
 # Read each name layer of the satellite data
 SRNames = st_layers(SR_GPKG)$name
 SR = lapply(SRNames, function(name) st_read(SR_GPKG, layer=name))
 # Remove geom column from the dataset
 SR_nogeom <- lapply(SR, function(x) sf::st_drop_geometry(x))

 # run a forloop through all the band layer

 for (layer in seq_along(SR)) {
    band_filled <- trans_multi_SR(SR_nogeom[[layer]]) # Interpolate missing data using stlplus

    df <- as.data.frame(band_filled$origin[[1]]) # extract colnames of the ts list
    df$date <- rownames(df)  # Add row names as a new column
    rownames(df) <- NULL     # Optional: Remove row names for cleaner data frame
    ts_notcomb <- df$date
    # create an empty list to collect interpolated sample ids
    tslist_SR_new <- list()
    for (row in 1:length(band_filled$filled)) {
      if(!row %in% remov_comb_ndvi) { tslist_SR_new <- append(tslist_SR_new,list(band_filled$filled[[row]]))}
    }
    # convert the final list into dataframe
    tslist_SR_new_df <- do.call(rbind.data.frame, tslist_SR_new)
    colnames(tslist_SR_new_df) <- df$date
    # Write the interpolated data layer with stlplus
    OutFile = paste0(main_dir, "Data_DRMAT/SR_cal_filled.gpkg")

    st_write(tslist_SR_new_df, dsn = OutFile, layer = SRNames[layer])
    # linear interpolation for SR cal. Combine stlplus + linear data
    SR_lininter <- linear_inter(tslist_non_and_lin_inter$remain_id, band_filled$origin, tslist_SR_new_df, ndvi = FALSE)

    # Write all interpolated data
    OutFile_SR = paste0(main_dir, "Data_DRMAT/SR_cal_filled_all.gpkg")
    st_write(SR_lininter, dsn = OutFile_SR, layer = SRNames[layer])
 }
 
 # save sample id that got deleted for calibration data due to Nan and not suitable for interpolation
 SR_nogeom[[1]] <- tibble::rowid_to_column(SR_nogeom[[1]], "ID")
 SR_nogeom[[1]] <- SR_nogeom[[1]][SR_nogeom[[1]]$ID %in% remov_comb_ndvi,]
 remove_sample_id <- SR_nogeom[[1]]$sample_id
 write.csv(remove_sample_id, paste0(main_dir,"./Data_DRMAT/deleted_sampleid_set.csv"))

# Prepare data for DRMAT format
 # read data filled by stlplus
 SR_GPKG_filled = "./Data/Data_DRMAT/SR_cal_filled.gpkg"
 SRNames = st_layers(SR_GPKG_filled)$name
 # read name layers
 SR_filled = lapply(SRNames, function(name) st_read(SR_GPKG_filled, layer=name))
# merge six bands and convert to a matrix of 3 dimensions
 SR_filled_3D <- abind(SR_filled[[2]], SR_filled[[3]], SR_filled[[4]], SR_filled[[5]], SR_filled[[6]], SR_filled[[7]], along=3)
# write the data
 write.csv(SR_filled_3D, paste0(main_dir, "Data_DRMAT/SR_3d_cal_filled.csv"))

 
 # read data filled by stlplus + linear interpolation 
 SR_GPKG_filled_all = "./Data/Data_DRMAT/SR_cal_filled_all.gpkg"
 # read name layers
 SR_filled_all = lapply(SRNames, function(name) st_read(SR_GPKG_filled_all, layer=name))
 # merge six bands and convert to a matrix of 3 dimensions
 SR_filled_3D_all <- abind(SR_filled_all[[2]], SR_filled_all[[3]], SR_filled_all[[4]], SR_filled_all[[5]], SR_filled_all[[6]], SR_filled_all[[7]], along=3)
 # write the data
 write.csv(SR_filled_3D_all, paste0(main_dir, "Data_DRMAT/SR_3d_cal_filled_all.csv"))
 ## Check whether the length of ndvi, sample_id and SR the same, similarly check the length of ymd and t
 if (nrow(tslist_non_and_lin_inter$tslist_all) == nrow(SR_filled_3D_all) == length(sample_id_ndvi_all_cal) && length(ymd_ndvi_cal) == length(t_ndvi_cal)){
   print("all files have expected length") else {
     print("one or more of the files does not have the expected length")
   }
 }

############################# USING STLPLUS TO INTERPOLATE MISSING DATA FOR NDVI VAL #############################
 # Read ndvi files for calibration dataset with cloudfree data
 l8ndvi_val <- st_read("./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_val.gpkg")
 # Interpolate missing data using stlplus
 ndvid_filled_val <- trans_multi_SR("./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_val.gpkg")
 # extract the results of the interpolated data: 1, a list of id that will not be interpolated and 2, ts with interpolated values
 remov_comb_ndvi_val <- ndvid_filled_val$rm_id
  tslist_filled_val <- ndvid_filled_val$filled

  # These steps help to track which sample id index is interpolated and which ones need to be removed
  # assign sample id again to the tslist that got interpolated
  # Track sample ids which ones are remained
  # turn ts lists into dataframe (result of interpolation)
 tslist_filled_df_val <- do.call(rbind.data.frame, tslist_filled_val)
 # bind sample id column to that ts df
 tslist_filled_df_val <- cbind(tslist_filled_df_val, sample_id = l8ndvi_val$sample_id)
 # create an id column for that ts df
 tslist_filled_df_val$id <- seq.int(nrow(tslist_filled_df_val))
 
 # remove the sample ids that could not be interpolated from the ts dataframe based in the id columns and rm_id output of the interpolation
 # create a ts df without all uninterpolated sample ids
 tslist_filled_df_rm_val <- tslist_filled_df_val[(!tslist_filled_df_val$id %in% remov_comb_ndvi_val),]
 
 
 # Append rows that were interpolated in a new list and later use this list for create final observation data
 for(row in 1:length(tslist_filled_val)) {if (length(tslist_filled_val[[row]]) == 186) {
   df <- as.data.frame(tslist_filled_val[[1]])
   df$date <- rownames(df)  # Add row names as a new column
   rownames(df) <- NULL     # Optional: Remove row names for cleaner data frame
   ts_notcomb <- df$date 
   break
 }}
 
 # assign filled ts list (all sample ids) to a new  variable 
 tslist_ndvi_186_val <- tslist_filled_val
 # create an empty list to collect interpolated sample ids
 tslist_new_val <- list()
 for (row in 1:length(tslist_ndvi_186_val)) {# only row index that not in the remove id list is selected
   # NOTE: this step seems to similar to the step done for the tslist_filled_df_rm. check
   if(!row %in% remov_comb_ndvi_val) { tslist_new_val <- append(tslist_new_val,list(tslist_ndvi_186_val[[row]]))}
 }
 # convert the final list into dataframe
 tslist_new_df_val <- do.call(rbind.data.frame, tslist_new_val)
 colnames(tslist_new_df_val) <- df$date
 # prepare sample id, time, ymd for DRMAT
 # Use tslist_filled_df which contain remaining sample ide for creating a sample id dataframe
 # created a dataframe of sample_ids that got interpolated with stlplus
 sample_id_ndvi_val <- get_sampleid(tslist_filled_df_rm_val, "Data_DRMAT/id_ndvi_val.csv")
 colnames(sample_id_ndvi_val) <- "sample_id"
 write.csv(sample_id_ndvi_val, paste0(main_dir, "Data_DRMAT/id_ndvi_val.csv"))
 # create ymd of the time series dates (186)
 ymd_ndvi_val <- getymd(tslist_new_df_val, "Data_DRMAT/ymd_ndvi_val.csv" )
 # create a dataframe with a fractional year for each date in the time series dates
 t_ndvi_val <- getfracyear(tslist_new_df_val, "Data_DRMAT/t_ndvi_val.csv")
 # write the final observation data with interpolated values done by stlplus
 outputname_ndvi_val <- paste0(main_dir, "Data_DRMAT/ndvi_val.csv")
 write.csv(tslist_new_df_val, outputname_ndvi_val)
 
 colnames(sample_id_ndvi_val) <- "sample_id"
ymd_test <- read.csv("/home/duyen/Master_Thesis/Data/Data_DRMAT/ymd_ndvi_val.csv", header = FALSE)
 ######################## USING LINEAR INTERPOLATION IN NDVI DATA ###############
 # interpolate data that was not interpolated by stlplus
 tslist_non_and_lin_inter_val <- linear_inter(remov_comb_ndvi_val, ndvid_filled_val$origin, tslist_new_df_val, ndvi = TRUE)

 ### adding id of linear interpolation
 # if (nrow(tslist_filled_df_val) == 19251) {
 #   tslist_filled_df_val$id <- seq.int(nrow(tslist_filled_df_val))} else {
 #     print("Total number of the observation is not equal 19251")
 #   }
 # get sample ids that were selected to be linearly interpolated
 sample_id_lininter_val <- as.data.frame((tslist_filled_df_val[tslist_filled_df_val$id %in% tslist_non_and_lin_inter_val$remain_id, ]$sample_id))
 colnames(sample_id_lininter_val) <- c("sample_id")
 # sample ids of stlplus
 colnames(sample_id_ndvi_val) <- c("sample_id")
 # append sample id of stlplus and linear interpolation
 sample_id_ndvi_all_val <- as.data.frame(c(as.integer(sample_id_ndvi_val$sample_id), sample_id_lininter_val$sample_id))
 colnames(sample_id_ndvi_all_val) <- c("sample_id")
 sample_id_ndvi_all_val <- as.integer(sample_id_ndvi_all_val$sample_id)
 # write new sample id for stlplus+linear (all-except for all nans?) and new ndvi cal (all interpolated-except for all nans?)
  write.csv(sample_id_ndvi_all_val, paste0(main_dir, "Data_DRMAT/id_ndvi_all_val.csv"))
 write.csv(tslist_non_and_lin_inter_val$tslist_all, paste0(main_dir, "Data_DRMAT/ndvi_all_val.csv"))
 
 ############################## USING STLPLUS + LINEAR TO INTERPOLATE MISSING DATA FOR SR VAL##########################
 # Read the validation set of satellite data 
 SR_GPKG_val = "./Intermediate product/cloud_free_product/__cloudfree_L8TS__SR_val.gpkg"
 # Read each name layer of the satellite data
 SRNames_val = st_layers(SR_GPKG_val)$name
 SR_val = lapply(SRNames_val, function(name) st_read(SR_GPKG_val, layer=name))
 # Remove geom column from the dataset
 SR_nogeom_val <- lapply(SR_val, function(x) sf::st_drop_geometry(x))
 
 # run a forloop through all the band layer
 check_list <- list()
 for (layer in seq_along(SR_val)) {
   band_filled_val <- trans_multi_SR(SR_nogeom_val[[layer]]) # Interpolate missing data using stlplus
   
   df <- as.data.frame(band_filled_val$origin[[1]])  # extract colnames of the ts list
   df$date <- rownames(df)  # Add row names as a new column
   rownames(df) <- NULL     # Optional: Remove row names for cleaner data frame
   ts_notcomb <- df$date
   # create an empty list to collect interpolated sample ids
   tslist_SR_new_val <- list()
   for (row in 1:length(band_filled_val$filled)) {
     if(!row %in% remov_comb_ndvi_val) { tslist_SR_new_val <- append(tslist_SR_new_val,list(band_filled_val$filled[[row]]))}
   }
   # convert the final list into dataframe
   tslist_SR_new_df_val <- do.call(rbind.data.frame, tslist_SR_new_val)
   colnames(tslist_SR_new_df_val) <- df$date
   # Write the interpolated data layer with stlplus
   OutFile_val = paste0(main_dir, "Data_DRMAT/SR_val_filled.gpkg")
   st_write(tslist_SR_new_df_val, dsn = OutFile_val, layer = SRNames_val[layer])
   # linear interpolation for SR cal. Combine stlplus + linear data
   SR_lininter_val <- linear_inter(tslist_non_and_lin_inter_val$remain_id, band_filled_val$origin, tslist_SR_new_df_val, ndvi = FALSE)
   check_list <- append(check_list, list(SR_lininter_val))
   # Write all interpolated data
   OutFile_SR_val = paste0(main_dir, "Data_DRMAT/SR_val_filled_all.gpkg")
   st_write(SR_lininter_val, dsn = OutFile_SR_val, layer = SRNames_val[layer])
   
   
 }
 
 # Prepare data for DRMAT format
 # read data filled by stlplus
 SR_GPKG_filled_val = "./Data/Data_DRMAT/SR_val_filled.gpkg"
 # read name layers
 SR_filled_val = lapply(SRNames_val, function(name) st_read(SR_GPKG_filled_val, layer=name))
 # merge six bands and convert to a matrix of 3 dimensions
 SR_filled_3D_val <- abind(SR_filled_val[[2]], SR_filled_val[[3]], SR_filled_val[[4]], SR_filled_val[[5]], SR_filled_val[[6]], SR_filled_val[[7]], along=3)
 # write the data
 write.csv(SR_filled_3D_val, paste0(main_dir, "Data_DRMAT/SR_3d_val_filled.csv"))
 
 
 # read data filled by stlplus + linear interpolation 
 SR_GPKG_filled_all_val = "./Data/Data_DRMAT/SR_val_filled_all.gpkg"
 # read name layers
 SR_filled_all_val = lapply(SRNames_val, function(name) st_read(SR_GPKG_filled_all_val, layer=name))
 # merge six bands and convert to a matrix of 3 dimensions
 SR_filled_3D_all_val <- abind(SR_filled_all_val[[2]], SR_filled_all_val[[3]], SR_filled_all_val[[4]], SR_filled_all_val[[5]], SR_filled_all_val[[6]], SR_filled_all_val[[7]], along=3)
 # write the data
 write.csv(SR_filled_3D_all_val, paste0(main_dir, "Data_DRMAT/SR_3d_val_filled_all.csv"))
 ## Check whether the length of ndvi, sample_id and SR the same, similarly check the length of ymd and t
 if (nrow(tslist_non_and_lin_inter_val$tslist_all) == nrow(SR_filled_3D_all_val) == length(sample_id_ndvi_all_val) && length(ymd_ndvi_val) == length(t_ndvi_val)){
   print("all files have expected length") else {
     print("one or more of the files does not have the expected length")
   }
 }
 ######################## PREPARE SR_CAL FOR COLD ###############################
# read SR_call filled by stlplus and has 3 d format
cold_data <- read.csv(paste0(main_dir, "Data_DRMAT/SR_3d_cal_filled.csv"))
# remove rows containing Nan
row_w_nan <- cold_data[!complete.cases(cold_data),]
# extract ID that got removed
rm_ID <- row_w_nan$X
rm_ID <- as.data.frame(rm_ID)
write_csv(rm_ID,paste0(main_dir, "Data_COLD/rm_id_w_nan_af_stlinter_SR_cal.csv"))
# remove rows containing Nan
# NOTE: THIS WAS NOT USED FOR DRMAT, but CAN BE USED FOR DRMAT
cold_data_nonan <- na.omit(cold_data)
write.csv(cold_data_nonan, paste0(main_dir, "Data_COLD/SR_cal_cold.csv"))

# prepare data for COLD
SR_GPKG_filled = "./Data/Data_DRMAT/SR_cal_filled.gpkg"
SRNames = st_layers(SR_GPKG_filled)$name
# read name layers
SR_filled = lapply(SRNames, function(name) st_read(SR_GPKG_filled, layer=name))
SR_filled_ID <- lapply(SR_filled, function(x) tibble::rowid_to_column(x, "ID"))
# remove rows with NAN
SR_rm_ID <- lapply(SR_filled_ID, function(x) x[-rm_ID,])
# Prepare data for COLD
mainDir_cold <- "/home/duyen/Master_Thesis/"
subDir <- "Data/Data_COLD"
ifelse(!dir.exists(file.path(mainDir_cold, subDir)), dir.create(file.path(mainDir_cold, subDir)), FALSE)
write.csv(SR_rm_ID[[2]], paste0(mainDir_cold, "Data/Data_COLD/Lansat_8_SRB2.csv"))
write.csv(SR_rm_ID[[3]], paste0(mainDir_cold, "Data/Data_COLD/Lansat_8_SRB3.csv"))
write.csv(SR_rm_ID[[4]], paste0(mainDir_cold, "Data/Data_COLD/Lansat_8_SRB4.csv"))
write.csv(SR_rm_ID[[5]], paste0(mainDir_cold, "Data/Data_COLD/Lansat_8_SRB5.csv"))
write.csv(SR_rm_ID[[6]], paste0(mainDir_cold, "Data/Data_COLD/Lansat_8_SRB6.csv"))
write.csv(SR_rm_ID[[7]], paste0(mainDir_cold, "Data/Data_COLD/Lansat_8_SRB7.csv"))


############### USING LAPPLY #######################
 # Run BFAST Lite on the dataset
 check_list_ap <- list()
SR_filled_all <- lapply(SR_nogeom_val, function(x) {
  band_filled_val_ap <- trans_multi_SR(SR_nogeom_val)
  
  df_ap <- as.data.frame(band_filled_val_ap$origin[[1]])
  df_ap$date <- rownames(df_ap)  # Add row names as a new column
  rownames(df_ap) <- NULL     # Optional: Remove row names for cleaner data frame
  ts_notcomb <- df_ap$date
  tslist_SR_new_val_ap <- list()
  for (row in 1:length(band_filled_val_ap$filled)) {
    if(!row %in% remov_comb_ndvi_val) { tslist_SR_new_val_ap <- append(tslist_SR_new_val_ap,list(band_filled_val_ap$filled[[row]]))}
  }
  tslist_SR_new_df_val_ap <- do.call(rbind.data.frame, tslist_SR_new_val_ap)
  colnames(tslist_SR_new_df_val_ap) <- df$date
  # Write the output of BFAST lite
  OutFile_val = paste0(main_dir, "Data_DRMAT/SR_val_filled_ap.gpkg")
  
  st_write(tslist_SR_new_df_val_ap, dsn = OutFile_val, layer = SRNames_val[layer])
  
  SR_lininter_val_ap <- linear_inter(tslist_non_and_lin_inter_val$remain_id, band_filled_val_ap$origin, tslist_SR_new_df_val_ap, ndvi = FALSE)
  check_list_ap <- append(check_list_ap, list(SR_lininter_val_ap))
  # Write the output of BFAST lite
  OutFile_SR_val_ap = paste0(main_dir, "Data_DRMAT/SR_val_filled_all_ap.gpkg")
  
  st_write(SR_lininter_val_ap, dsn = OutFile_SR_val_ap, layer = SRNames_val[layer])
  my_file <- c( "nonlinfilled" = band_filled_val_ap, "linandnonlinfilled" = SR_lininter_val_ap )
  return(my_file)
})



length()

################ NAN stats ####################
# create a table that contains different consecutive nans scenario
list_frq <- c()
 for (row in  1:length(ndvid_filled_cal$origin)) {
   x = is.na(ndvid_filled_cal$origin[[row]])
   
   if (all(x) != TRUE) {
   x_1 <- rle(x)$lengths[rle(x)$values]
   X1_frq <- as.data.frame(table(x_1))
   list_frq <- c(list_frq, as.integer(as.character(X1_frq$x_1)))} else {
     list_frq <- c(list_frq, length(x))
   }
 }
frq_tb <- as.data.frame(table(list_frq))

# Barplot-plot all scenario's of consecutive nans numbers
ggplot(frq_tb, aes(x=list_frq, y=Freq)) + geom_bar(stat = "identity")

ggplot(frq_tb, aes(x=as.factor(list_frq), fill=as.factor(list_frq) )) + 
  geom_bar( ) +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position="none")



