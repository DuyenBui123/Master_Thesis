# ----  General info --------------------------
#' 
#' Name:             Duyen Bui
#' Study programs:   Master of Geo-information Science 
#' Student number:   1044767
#' 
#' Project:          Thesis: Compare the performance of multivariate and univariate time series models in detecting global land cover change
#' Supervisors:      Dainius Masiliunas & Nandin-Erdene Tsendbazar
#'
#' Script to preprocess time series data including: merge IIASA and WUR data, change the names of time series columns in all layers, and calculate NDVI index
#'
#'_____________________________________________________________________

#### Download and install packages ####
# LIBRARY
# create notin operator
`%notin%` <- Negate(`%in%`)

# Download packages if not available
pckgs <- c("tidyverse","sf", "zoo", "parallel","reshape2", "pbapply" )

if (any(pckgs %notin% rownames(installed.packages())==TRUE)){
  install.packages(pckgs, repos = c(CRAN = "http://cloud.r-project.org"))}

# LOAD PACKAGE
sapply(pckgs, FUN = require, character.only = TRUE)
library(sf)
library(zoo)
library(parallel)
library(reshape2)
library(pbapply)
library(tidyverse)

source("C:\\Master_Thesis\\GitHub_code\\postprocessing-bfast\\src\\utils\\utils.r")
source("C:\\Master_Thesis\\GitHub_code\\postprocessing-bfast\\src\\utils\\covariate-names.r")
source("C:\\Master_Thesis\\GitHub_code\\postprocessing-bfast\\src\\utils\\load-sampling-data.r")
source("C:\\Master_Thesis\\Code\\Utils\\preprocessing_ref_and_sat_data.R")

# Assign input path to variable

SR_GPKG_WUR = "C:\\Master_Thesis\\Data\\WURChange20152019_Landsat8_TS.gpkg"
SR_GPKG_IIASA = "C:\\Master_Thesis\\Data\\IIASAChange20152018_Landsat8_TS.gpkg"

# Assign output path to variable
Feature_GPKG_val = "C:\\Master_Thesis\\Data\\validation_features\\"
Feature_GPKG_cal = "C:\\Master_Thesis\\Data\\calibration_features\\"

# Read csv file of calibration and validation data
cal_compl_data <- read.csv("C:\\Master_Thesis\\Intermediate product\\cal_compl_data.csv")
val_compl_data <- read.csv("C:\\Master_Thesis\\Intermediate product\\val_compl_data.csv")
################################# MERGE WUR AND IIASA SATALLITE DATA #######################################################
# Read all band layer of satelite WUR and IIASA datalayers into a single list "SR"
SRNames_WUR = st_layers(SR_GPKG_WUR)$name
SR_WUR = lapply(SRNames_WUR, function(name) st_read(SR_GPKG_WUR, layer=name)) # SR is a flatten satelite collection, each SR is a band and each band contains all the sample ids

SRNames_IIASA = st_layers(SR_GPKG_IIASA)$name
SR_IIASA <- lapply(SRNames_IIASA, function(name) st_read(SR_GPKG_IIASA, layer=name)) # SR is a flatten satelite collection, each SR is a band and each band contains all the sample ids

# Merge two satelite dataset for validation data together and align the satelite data with the reference data to remove sample ids that are not included in reference data
SR_val <- merge_sat_data(SR_WUR, SR_IIASA, val_compl_data)

# Write a CSV for the merged validation data
OutFile = paste(Feature_GPKG_val, "SR_val.gpkg", sep="_")


if (!file.exists(OutFile)) {
  for (layer in seq_along(SR_val)) { 
    st_write(SR_val[[layer]], dsn = OutFile, layer = SRNames_IIASA[layer])}}
# Merge two satelite dataset for calibration data together and align the satelite data with the reference data to remove sample ids that are not included in reference data
SR_cal <- merge_sat_data(SR_WUR, SR_IIASA, cal_compl_data)

# Write a CSV for the merged calibration data
OutFile_cal = paste(Feature_GPKG_cal, "SR_cal.gpkg", sep="_")

if (!file.exists(OutFile_cal)) {
  for (layer in seq_along(SR_cal)) { 
    st_write(SR_cal[[layer]], dsn = OutFile_cal, layer = SRNames_IIASA[layer])}}

# Check whether the sample_id of the merged calibration and merged validation are dublicated or not
lapply(SR_val, function(x, y) {if (any(x$sample_id %in% y$sample_id)) {print("There are dublicated between cal and val data")} else {print("no dublicated found")}}, SR_cal)

################################# RENAME TIME SERIES COLUMND, CALCULATE NDVI #######################################################
# For calibration data
# Save a template for writing, with consistent column names for calibration data
OutTemplate_cal = SR_cal[[1]]
names(OutTemplate_cal)[datecols(OutTemplate_cal)] = strtrim(names(OutTemplate_cal)[datecols(OutTemplate_cal)], 11)
# Extract the time series in matrix format from an SF object by date
# Convert into a zoo matrix to correct the time series column names
SRZ_cal = lapply(lapply(SR_cal, SFToMatrix), MatrixToZoo)
names(SRZ_cal) = SRNames_WUR
# Calculate NDVI index
NDVI  = function(BLUE=NULL, GREEN=NULL, RED     , NIR,  SWIR=NULL)
{return((NIR-RED)/(NIR+RED))}
# Write a csv for NDVI calibration data
OutFile = paste(Feature_GPKG_cal, "NDVI_cal.gpkg", sep="_")
if (!file.exists(OutFile)) {
  NDVIz = NDVI(SRZ_cal[["SR_B2"]], SRZ_cal[["SR_B3"]], SRZ_cal[["SR_B4"]], SRZ_cal[["SR_B5"]], SRZ_cal[["SR_B6"]])
  NDVIsf = ZooToSF(NDVIz, OutTemplate_cal)
  print(NDVIz)
  st_write(NDVIsf, OutFile)
  rm(NDVIsf, NDVIz)
}

# For validation data
# Save a template for writing, with consistent column names for validation data
OutTemplate_val = SR_val[[1]]
names(OutTemplate_val)[datecols(OutTemplate_val)] = strtrim(names(OutTemplate_val)[datecols(OutTemplate_val)], 11)

# Convert into a zoo matrix to correct the time series column names
SRZ_val = lapply(lapply(SR_val, SFToMatrix), MatrixToZoo)
names(SRZ_val) = SRNames_WUR
# Calculate all indices of interest
NDVI  = function(BLUE=NULL, GREEN=NULL, RED     , NIR,  SWIR=NULL)
{return((NIR-RED)/(NIR+RED))}
# Write a csv for NDVI validation data
OutFile = paste(Feature_GPKG_val, "NDVI_val.gpkg", sep="_")
if (!file.exists(OutFile)) {
  NDVIz = NDVI(SRZ_val[["SR_B2"]], SRZ_val[["SR_B3"]], SRZ_val[["SR_B4"]], SRZ_val[["SR_B5"]], SRZ_val[["SR_B6"]])
  NDVIsf = ZooToSF(NDVIz, OutTemplate_val)
  st_write(NDVIsf, OutFile)
  rm(NDVIsf, NDVIz)
}

# Remove the surfix (_SR_...) in the time series column names for all the layers in the satellite data, so it only contain year, month, day with this format X????.??.??
# For calibration data
for (layer in seq_along(SR_cal)) {
  SR_cal[[layer]] <- ZooToSF(SRZ_cal[[layer]], OutTemplate_cal)
}
# Write CSV for the calibration satellite data
OutFile_cal_cor = paste(Feature_GPKG_cal, "_cor_SR_cal.gpkg", sep="_")

if (!file.exists(OutFile_cal_cor)) {
  for (layer in seq_along(SR_cal)) { 
    st_write(SR_cal[[layer]], dsn = OutFile_cal_cor, layer = SRNames_IIASA[layer])}}

# For validation data

for (layer in seq_along(SR_val)) {
  SR_val[[layer]] <- ZooToSF(SRZ_val[[layer]], OutTemplate_val)
}

# Write CSV for the validation satellite data
OutFile_val_cor = paste(Feature_GPKG_val, "_cor_SR_val.gpkg", sep="_")

if (!file.exists(OutFile_val_cor)) {
  for (layer in seq_along(SR_val)) { 
    st_write(SR_val[[layer]], dsn = OutFile_val_cor, layer = SRNames_IIASA[layer])}}

             