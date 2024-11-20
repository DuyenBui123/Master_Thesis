
# Load library 
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(abind, stlplus, TSstudio, zoo, here, devtools, sf, data.table, tidyverse, pbapply, parallel, data.table)
install.packages("imputeTS")
library(imputeTS)
source("/home/duyen/Master_Thesis/GitHub_code/BFASTutils.R")
source("/home/duyen/Master_Thesis/GitHub_code/getTrimmedts.R")
source("/home/duyen/Master_Thesis/Code/Utils/Utils.R")
setwd("/home/duyen/Master_Thesis/")
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

# # Read the CALIBRATION set of satellite data 
# SR_GPKG = "./Intermediate product/cloud_free_product/__cloudfree_L8TS__SR_cal.gpkg"
# # Read each name layer of the satellite data
# SRNames = st_layers(SR_GPKG)$name
# SR = lapply(SRNames, function(name) st_read(SR_GPKG, layer=name))
# 
# 
# # write dates
# outputname_ymd <- paste0(main_dir, "Data_DRMAT/ymd.csv")
# date_ts <- colnames(SR_df[[1]])
# year_ts <- c()
# for (index in 1:length(date_ts)) { year_ts <- append(year_ts,substring(date_ts[index], 2, 5))
# }
# year_ts_df <- as.data.frame(year_ts)
# 
# 
# month_ts <- c()
# for (index in 1:length(date_ts)) { month_ts <- append(month_ts,substring(date_ts[index], 7, 8))
# }
# month_ts_df <- as.data.frame(month_ts)
# 
# 
# day_ts <- c()
# for (index in 1:length(day_ts)) { day_ts <- append(day_ts,substring(date_ts[index], 10, 11))
# }
# day_ts_df <- as.data.frame(day_ts)
# 
# date_ts_df <- cbind(year_ts_df, month_ts_df, day_ts_df)
# colnames(date_ts_df) <- NULL
# write.csv(date_ts_df, outputname_ymd)
# # Write id
# SR_df_id <- lapply(SR,function(x) x <- x[, !names(x) %in% c("tile","centroid_x", "centroid_y")] %>% st_drop_geometry())
# outputname <- paste0(main_dir, "Data_DRMAT/id.csv")
# sample_id <- as.integer(SR_df_id[[1]][, names(SR_df_id[[1]]) %in% c("sample_id")])
# sample_id <- as.data.frame(sample_id)
# colnames(sample_id) <- NULL
# write.csv(sample_id, outputname)
# 
# # Write SR
# 
# OutFile = paste0(main_dir, "Data_DRMAT/SR_3D.Rdata")
# OutFile_rds = paste0(main_dir, "Data_DRMAT/SR_3D.rds")
# OutFile_mat = paste0(main_dir, "Data_DRMAT/SR_3D.mat" )
# SR_df <- lapply(SR,function(x) x <- x[, !names(x) %in% c("tile", "sample_id","centroid_x", "centroid_y")] %>% st_drop_geometry())
# 
# SR_3D <- abind(SR_df[[2]], SR_df[[3]], SR_df[[4]], SR_df[[5]], SR_df[[6]], SR_df[[7]], along=3)
# SR_3D <- unname(SR_3D)
# dim(SR_3D)
# # Interpolate
# trans_SR_3D <- lapply(SR_3D, function (x) t(x))
# 
# ###################
# # write fractional year
# date_ts <- sub('X','',date_ts)
# date_ts <-  gsub("\\.", "-", date_ts, perl = TRUE)
# t <- c()
# for (index in 1:length(date_ts)) {
#   t <- append(t, sprintf("%.11f",decimal_date(as.POSIXlt(date_ts[index]))))
# }
# options(digits=11)
# t <- as.numeric(t)
# t <- as.data.frame(t)
# colnames(t) <- NULL
# outputname_t <- paste0(main_dir, "/Data_DRMAT/t.csv")
# write.csv(t, outputname_t)
###################### prepare for ndvi_cal #################################
# write nivd
pbo <- pboptions(type="timer")
mycores <- detectCores()
ndvid_filled_cal <- trans_multi_SR("./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg")
remov_comb <- append(remov_index, diff_length)
tslist_filled <- ndvid_filled_cal$filled
# Plot time series
test <- ts(tslist_spid2, start  = c(2013, 4), end = c(2021, 12), frequency = 22)
ts.plot(test,tslist[[2]],
        col = 2:1)

# prepare the NDVI for DRMAT
# assign sample id again to the tslist got interpolated

# Track sample id which one is remained
tslist_filled_df <- do.call(rbind.data.frame, tslist_filled)
tslist_filled_df <- cbind(tslist_filled_df, sample_id = l8ndvi$sample_id)
tslist_filled_df$id <- seq.int(nrow(tslist_filled_df))
# remove the sample id that could not be interpolated
tslist_filled_df_rm <- tslist_filled_df[(!tslist_filled_df$id %in% c(remov_index, diff_length)),]


# Append rows that were interpolated in a new list and later use this list for create final observation data
for(row in 1:length(tslist_filled)) {if (length(tslist_filled[[row]]) == 186) {
  df <- as.data.frame(tslist_filled[[1]])
  df$date <- rownames(df)  # Add row names as a new column
  rownames(df) <- NULL     # Optional: Remove row names for cleaner data frame
  ts_notcomb <- df$date 
  break
}}


tslist_ndvi_186 <- tslist_filled
tslist_new <- list()
for (row in 1:length(tslist_ndvi_186)) {
  if(!row %in% remov_comb) { tslist_new <- append(tslist_new,list(tslist_ndvi_186[[row]]))}
}
tslist_new_df <- do.call(rbind.data.frame, tslist_new)
colnames(tslist_new_df) <- df$date
# prepare sample id, time, ymd for DRMAT

# Use tslist_filled_df which contain remaining sample ide for creating a sample id dataframe
sample_id_ndvi_cal <- get_sampleid(tslist_filled_df_rm, "Data_DRMAT/id_ndvi_cal.csv")
ymd_ndvi_cal <- getymd(tslist_new_df, "Data_DRMAT/ymd_ndvi_cal.csv" )
t_ndvi_cal <- getfracyear(tslist_new_df, "Data_DRMAT/t_ndvi_cal.csv")

outputname_ndvi_cal <- paste0(main_dir, "Data_DRMAT/ndvi_cal.csv")
write.csv(tslist_new_df, outputname_ndvi_cal)

############################## PREPARE DATA FOR SR CAL##########################
# Read the CALIBRATION set of satellite data 
SR_GPKG = "./Intermediate product/cloud_free_product/__cloudfree_L8TS__SR_cal.gpkg"
# Read each name layer of the satellite data
SRNames = st_layers(SR_GPKG)$name
SR = lapply(SRNames, function(name) st_read(SR_GPKG, layer=name))
# Remove geom column from the dataset
SR_nogeom <- lapply(SR, function(x) sf::st_drop_geometry(x))

# Run BFAST Lite on the dataset
for (layer in seq_along(SR)) {
  band_filled <- trans_multi_SR(SR_nogeom[[layer]])
  df <- as.data.frame(band_filled$origin[[1]])
  df$date <- rownames(df)  # Add row names as a new column
  rownames(df) <- NULL     # Optional: Remove row names for cleaner data frame
  ts_notcomb <- df$date
  tslist_new <- list()
  for (row in 1:length(band_filled$filled)) {
    if(!row %in% remov_comb) { tslist_new <- append(tslist_new,list(band_filled$filled[[row]]))}
  }
  tslist_new_df <- do.call(rbind.data.frame, tslist_new)
  colnames(tslist_new_df) <- df$date
  # Write the output of BFAST lite
  OutFile = paste0(main_dir, "Data_DRMAT/SR_cal_filled.gpkg")
  
  st_write(tslist_new_df, dsn = OutFile, layer = SRNames[layer])
  
}



SR_GPKG_filled = "./Data/Data_DRMAT/SR_cal_filled.gpkg"
OutFile = paste0(main_dir, "Data_DRMAT/SR_3D.Rdata")
SR_filled = lapply(SRNames, function(name) st_read(SR_GPKG_filled, layer=name))

library(abind)
SR_filled_3D <- abind(SR_filled[[2]], SR_filled[[3]], SR_filled[[4]], SR_filled[[5]], SR_filled[[6]], SR_filled[[7]], along=3)
SR_3D <- unname(SR_3D)
write.csv(SR_filled_3D, paste0(main_dir, "Data_DRMAT/SR_3d_cal_filled.csv"))

test <- read.csv(paste0(main_dir, "Data_DRMAT/SR_3d_cal_filled.csv"))



############### USING LAPPLY #######################
# Run BFAST Lite on the dataset
SR_wonan <- lapply(SR_nogeom[[layer]], function(x) {
  band_filled <- trans_multi_SR(SR_nogeom[[layer]])
  df <- as.data.frame(band_filled$origin[[1]])
  df$date <- rownames(df)  # Add row names as a new column
  rownames(df) <- NULL     # Optional: Remove row names for cleaner data frame
  ts_notcomb <- df$date
  tslist_new <- list()
  for (row in 1:length(band_filled$filled)) {
    if(!row %in% remov_comb) { tslist_new <- append(tslist_new,list(band_filled$filled[[row]]))}
  }
  tslist_new_df <- do.call(rbind.data.frame, tslist_new)
  colnames(tslist_new_df) <- df$date
  return(tslist_new_df)
})

SR_nogeom[[1]] <- tibble::rowid_to_column(SR_nogeom[[1]], "ID")
SR_nogeom[[1]] <- SR_nogeom[[1]][SR_nogeom[[1]]$ID %in% remov_comb,]
remove_sample_id <- SR_nogeom[[1]]$sample_id
write.csv(remove_sample_id, paste0(main_dir,"./Data_DRMAT/deleted_sampleid_set.csv"))



################ NAN stats ####################
plot(which(is.na(ndvid_filled_cal$origin[[1]])), xlim = c(2013,2020), ylim = c(0,200), 
     ylab = "Index of missing value within dataset", xlab = "Index within the missing value")
ndvid_filled_cal$origin[[1]]










# set.seed(42)  # For reproducibility
# length(tslist[[5]])
# ts_select_1 <- window(tslist[[3]], start = 2014, end = 2020, frequency = 22)
# length(ts_select_1)
# # Parameters
# n <- 133      # Total number of observations
# n.p <- 22     # Seasonal period (e.g., weekly data in a year)
# 
# # Generate a vector with a high proportion of NA
# Y <- sample(c(NA, rnorm(1, mean = 50, sd = 10)), size = n, replace = TRUE, prob = c(0.85, 0.15))
# Y
# 
# # Generate cycle indices for the seasonal period
# cycleSubIndices <- rep(1:n.p, ceiling(n / n.p))[1:n]
# 
# # Check for fully missing subseries
# if (any(by(ts_select_1, list(cycleSubIndices), function(x) all(is.na(x))))) {
#   stop("There is at least one subseries for which all values are missing.")
# }
# length(ts_select_1)
# # Print results
# print("Generated Y:")
# print(Y)
# 
# print("Cycle Sub Indices:")
# print(cycleSubIndices)
# 
# print("Subseries grouped by cycle indices:")
# print(by(ts_select_1, list(cycleSubIndices), function(x) x))
# 
# # Define cycleSubIndices
# cycleSubIndices <- rep(1:n.p, ceiling(length(ts_select_1) / n.p))[1:length(ts_select_1)]
# 
# # Identify subseries with all NA values and print the subseries values
# missing_subseries <- by(ts_select_1, list(cycleSubIndices), function(x) {
#   
#   # Print the subseries being inspected
#   cat("Inspecting subseries:\n")
#   print(x)
#   
#   # Check if all values are NA
#   return(all(is.na(x)))
# })
# 
# # Check which subseries have all NA values and print their indices
# subseries_with_na <- which(missing_subseries)  # Returns the indices of subseries with all NA values
# 
# # If there are subseries with all NA, print their indices
# if (length(subseries_with_na) > 0) {
#   cat("The following subseries contain all missing data (NA):\n")
#   print(subseries_with_na)
# } else {
#   cat("No subseries contain all missing data (NA).")
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
