# ----  General info --------------------------
#' 
#' Name:             Duyen Bui
#' Study programs:   Master of Geo-information Science 
#' Student number:   1044767
#' 
#' Project:          Thesis: Compare the performance of multivariate and univariate time series models in detecting global land cover change
#' Supervisors:      Dainius Masiliunas & Nandin-Erdene Tsendbazar
#'
#' Script to filter out cloud for NDVI, and satelite data
#'
#'_____________________________________________________________________
# Load libraries -------------------------------------------------
# load packages

# Download packages if not available
options(repos = getOption("repos")["CRAN"])
if (!require("pacman")) install.packages("pacman"); library(pacman)
pckgs <- c("here", "sf", "zoo", "devtools","tidyverse", "pbapply", "lubridate" )
install.packages("devtools")
if (any(pckgs %notin% rownames(installed.packages())==TRUE)){
  install.packages(pckgs, repos = c(CRAN = "http://cloud.r-project.org"))}

install.packages("xtable")
sapply(pckgs, FUN = require, character.only = TRUE)

#source external functions
here::i_am("Code\\filter_time_series.R")
source(here("GitHub_code","postprocessing-bfast", "src","utils", "utils.r"))
source(here("GitHub_code","postprocessing-bfast", "src","utils", "covariate-names.r"))
source(here("GitHub_code", "postprocessing-bfast", "src", "utils", "load-sampling-data.r"))
source(here("Code","Utils", "Utils.r"))

# input path
SR_GPKG_CAL = "C:/Master_Thesis/Data/calibration_features/__cor_SR_cal.gpkg"
SR_GPKG_VAL = "C:/Master_Thesis/Data/validation_features/__cor_SR_val.gpkg"

# Output path
NDVI_cal_path <- ("C:/Master_Thesis/Data/calibration_features/_NDVI_cal.gpkg")
NDVI_val_path <- ("C:/Master_Thesis/Data/validation_features/_NDVI_val.gpkg")
# determine the name of the dataset 
datasetname <- ("C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\")

devtools::source_url("https://raw.githubusercontent.com/GreatEmerald/master-classification/master/src/pixel-based/utils/ProbaVDataDirs.r")
devtools::source_url("https://raw.githubusercontent.com/GreatEmerald/master-classification/master/src/pixel-based/utils/db-io.r")

############################################################################################
# Cloud removal for NDVI calibration dataset using LOESS function 
cor_SRNames_cal = st_layers(SR_GPKG_CAL)$name
cor_SR_cal <-  lapply(cor_SRNames_cal, function(name) st_read(SR_GPKG_CAL, layer=name)) # Read the corrected satellite calibration dataset
cloudfreeNDVI_cal <- filter_time_series(cor_SR_cal, SR_GPKG_CAL, NDVI_cal_path) # Cloud removal
# Write the filtered NDVI observations to a data file
outfilename <- paste0(datasetname, "_cloudfree_L8TS_",
                      "NDVI_cal", ".gpkg")
st_write(cloudfreeNDVI_cal, outfilename, "NDVI_cal", driver = 'GPKG', update = TRUE)


# Cloud removal for NDVI validation dataset using LOESS function
cor_SR_val <-  lapply(cor_SRNames_cal, function(name) st_read(SR_GPKG_VAL, layer=name))# Read the corrected satellite validation dataset
cloudfreeNDVI_val <- filter_time_series(cor_SR_val, SR_GPKG_VAL, NDVI_val_path) # Cloud removal
# Write the filtered NDVI validation to a data file 
outfilename <- paste0(datasetname, "_cloudfree_L8TS_",
                      "NDVI_val", ".gpkg")
st_write(cloudfreeNDVI_val, outfilename, driver = 'GPKG')

# For calibration satellite data
cor_SRNames_cal = st_layers(SR_GPKG_CAL)$name
cloudfree_SR_cal <- lapply(cor_SR_cal, filter_time_series, SR_GPKG_CAL, "none") # Cloud removal
# Write the filtered calibration satellite data to a data file 
OutFile_cal = paste(datasetname, "_cloudfree_L8TS_", "SR_cal.gpkg", sep="_")
if (!file.exists(OutFile_cal)) {
  for (layer in seq_along(cor_SR_cal)) { 
    st_write(cloudfree_SR_cal[[layer]], dsn = OutFile_cal, layer = cor_SRNames_cal[layer])}}


# For validation satellite data
cloudfree_SR_val <- lapply(cor_SR_val, filter_time_series, SR_GPKG_VAL, "none")
# Write the filtered calibration satellite data to a data file 
OutFile_val = paste(datasetname, "_cloudfree_L8TS_", "SR_val.gpkg", sep="_")
if (!file.exists(OutFile_val)) {
  for (layer in seq_along(cor_SR_cal)) { 
    st_write(cloudfree_SR_val[[layer]], dsn = OutFile_val, layer = cor_SRNames_cal[layer])}}

