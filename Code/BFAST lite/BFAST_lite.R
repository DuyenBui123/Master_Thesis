# ----  General info --------------------------
#' 
#' Name:             Duyen Bui
#' Study programs:   Master of Geo-information Science 
#' Student number:   1044767
#' 
#' Project:          Thesis: Compare the performance of multivariate and univariate time series models in detecting global land cover change
#' Supervisors:      Dainius Masiliunas & Nandin-Erdene Tsendbazar
#'
#' Script to calibrate BFAST lite on NDVI, and satellite data. 
#'
#'_____________________________________________________________________

install.packages("pacman")
library(pacman)
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(here, sf, zoo, data.table, tidyverse, pbapply,parallel, sandwich) #, bfast, strucchangeRcpp)

# source external functions
setwd("/home/duyen/Master_Thesis/")
# setwd("/home/duyen/Master_Thesis/GitHub_code/strucchangeRcpp/")
# setwd("/home/duyen/Master_Thesis/GitHub_code/bfast/")
here::i_am("Code/BFAST_lite.R")
source(here("GitHub_code", "Utils.R"))
source(here::here("GitHub_code", "getTrimmedts.R"))
debugSource(here::here("GitHub_code", "BFASTutils.R"))
source(here::here("GitHub_code", "cglops-change-detection", "src", "bfast-cal", "04-validation.r"))

# add progress bar option to show it in the terminal 

pbo <- pboptions(type="timer")
mycores <- detectCores()
set.seed(123, kind = "L'Ecuyer-CMRG" )
# set pathes and create folder
base_path <- "./Intermediate product/cloud_free_product/"

# Create output folder
mainDir <- "C:/Master_Thesis/"
subDir <- "Output"
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)
# load data ----------------------------------------------


# Calibrate BIC, formula: response ~ trend + harmon OR response ~ trend OR response ~ harmon, magnitude threshold: 0.25, 0.3, 1, -inf

BFASTlite_BIC_SandT_025 <- cal_BFAST("BIC", 0.25, 
                                       response ~ harmon + trend, "trend", "trend", 
                                       "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg", 
                                       "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_SandT_025.csv", plot = FALSE, cl=mycores)


BFASTlite_BIC_T_025 <- cal_BFAST("BIC", 0.25, response ~ trend, "trend", "trend",
                                   "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg", 
                                   "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_T_025.csv", plot = FALSE, cl=mycores)


BFASTlite_BIC_S_025 <- cal_BFAST("BIC", 0.25, response ~ harmon,
                                   c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                   c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                   "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg",
                                   "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_S_025.csv", plot = FALSE, cl=mycores)



BFASTlite_BIC_SandT_030 <- cal_BFAST("BIC", 0.30, 
                                     response ~ harmon + trend, "trend", "trend", 
                                     "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg", 
                                     "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_SandT_030.csv", plot = FALSE, cl=mycores)

BFASTlite_BIC_T_030 <- cal_BFAST("BIC", 0.30, response ~ trend, "trend", "trend",
                                 "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg", 
                                 "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_T_030.csv", plot = FALSE, cl=mycores)


BFASTlite_BIC_S_030 <- cal_BFAST("BIC", 0.30, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg",
                                 "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_S_030.csv", plot = FALSE, cl=mycores)

BFASTlite_BIC_SandT_050 <- cal_BFAST("BIC", 0.50, 
                                     response ~ harmon + trend, "trend", "trend", 
                                     "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg", 
                                     "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_SandT_050.csv", plot = FALSE, cl=mycores)

BFASTlite_BIC_T_050 <- cal_BFAST("BIC", 0.50, response ~ trend, "trend", "trend",
                                 "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg", 
                                 "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_T_050.csv", plot = FALSE, cl=mycores)

BFASTlite_BIC_S_050 <- cal_BFAST("BIC", 0.50, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg",
                                 "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_S_050.csv", plot = FALSE, cl=mycores)


BFASTlite_BIC_SandT_1 <- cal_BFAST("BIC", 1, response ~ harmon + trend, "trend", "trend",
                                   "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg", 
                                   "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_SandT_1.csv", plot = FALSE, cl=mycores)


BFASTlite_BIC_T_1 <- cal_BFAST("BIC", 1, response ~ trend, "trend", "trend",
                                 "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg", 
                                 "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_T_1.csv", plot = FALSE, cl = mycores)

BFASTlite_BIC_S_1 <- cal_BFAST("BIC", 1, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg",
                                 "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_S_1.csv", plot = FALSE, cl = mycores)



BFASTlite_BIC_SandT_inf <- cal_BFAST("BIC", -Inf, response ~ harmon + trend, "trend", "trend",
                                   "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg", 
                                   "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_SandT_Inf.csv", plot = FALSE, cl = mycores)

BFASTlite_BIC_T_inf <- cal_BFAST("BIC", -Inf, response ~ trend, "trend", "trend",
                               "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg", 
                               "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_T_Inf.csv", plot = FALSE, cl = mycores)


BFASTlite_BIC_S_inf <- cal_BFAST("BIC", -Inf, response ~ harmon,
                               c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                               c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                               "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg",
                               "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_S_Inf.csv", plot = FALSE, cl = mycores)

# Calibrate LWZ, formula: response ~ trend + harmon OR response ~ trend OR response ~ harmon, magnitude threshold: 0.25, 0.3, 1, -inf

BFASTlite_LWZ_SandT_025 <- cal_BFAST("LWZ", 0.25, 
                                     response ~ harmon + trend, "trend", "trend", 
                                     "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg", 
                                     "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_SandT_025.csv", plot = FALSE, cl = mycores)


BFASTlite_LWZ_T_025 <- cal_BFAST("LWZ", 0.25, response ~ trend, "trend", "trend",
                                 "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg", 
                                 "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_T_025.csv",plot = FALSE, cl = mycores)

BFASTlite_LWZ_S_025 <- cal_BFAST("LWZ", 0.25, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg",
                                 "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_S_025.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_SandT_030 <- cal_BFAST("LWZ", 0.30, 
                                     response ~ harmon + trend, "trend", "trend", 
                                     "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg", 
                                     "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_SandT_030.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_T_030 <- cal_BFAST("LWZ", 0.30, response ~ trend, "trend", "trend",
                                 "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg", 
                                 "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_T_030.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_S_030 <- cal_BFAST("LWZ", 0.30, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg",
                                 "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_S_030.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_SandT_050 <- cal_BFAST("LWZ", 0.50, 
                                     response ~ harmon + trend, "trend", "trend", 
                                     "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg", 
                                     "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_SandT_050.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_T_050 <- cal_BFAST("LWZ", 0.50, response ~ trend, "trend", "trend",
                                 "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg", 
                                 "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_T_050.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_S_050 <- cal_BFAST("LWZ", 0.50, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg",
                                 "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_S_050.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_SandT_1 <- cal_BFAST("LWZ", 1, response ~ harmon + trend, "trend", "trend", 
                                   "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg", 
                                   "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_SandT_1.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_T_1 <- cal_BFAST("LWZ", 1, response ~ trend, "trend", "trend",
                               "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg", 
                               "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_T_1.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_S_1 <- cal_BFAST("LWZ", 1, response ~ harmon,
                               c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                               c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                               "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg",
                               "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_S_1.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_SandT_inf <- cal_BFAST("LWZ", -Inf, response ~ harmon + trend, "trend", "trend", 
                                     "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg", 
                                     "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_SandT_inf.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_T_inf <- cal_BFAST("LWZ", -Inf, response ~ trend, "trend", "trend",
                                 "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg", 
                                 "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_T_inf.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_S_inf <- cal_BFAST("LWZ", -Inf, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg",
                                 "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_S_inf.csv", plot = FALSE, cl = mycores)



###################### Run BFAST lite ON NDVI VALIDATION DATASET ###############

# Run BFAST Lite on NDVI validation dataset
# Set of parameters is chosen based on F1 score of all calibration sets above. The one that has the highest F1 score is chosen
# BIC, response ~ trend, and magnitude threshold = 0.25 
BFASTlite_BIC_T_025_val <- cal_BFAST("BIC", 0.25, response ~ trend, "trend", "trend",
                                 "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_val.gpkg", 
                                 "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_T_025_val.csv",  plot = FALSE, cl = mycores)

####################### METHOD 1 for running multiple bands of satellite data ##
# Select the most dominant floor year that contain more than or equal to 3 breakpoints 
# from all bands. Then choose the highest magnitude breakpoints.
###################### RUN BFAST LITE ON SATELLITE DATA  #######################
# Read the CALIBRATION set of satellite data 
SR_GPKG = "./Intermediate product/cloud_free_product/__cloudfree_L8TS__SR_cal.gpkg"
# Read each name layer of the satellite data
SRNames = st_layers(SR_GPKG)$name
SR = lapply(SRNames, function(name) st_read(SR_GPKG, layer=name))
# Remove geom column from the dataset
SR_nogeom <- lapply(SR, function(x) sf::st_drop_geometry(x))
# Run BFAST Lite on the dataset
for (layer in seq_along(SR)) { breakpoint_SR <- cal_BFAST("BIC", 0.25, response ~ trend, "trend", "trend",
                                                          SR_nogeom[[layer]], 
                                         NULL,  plot = FALSE, cl = mycores)
                              # Write the output of BFAST lite
                              OutFile = paste(base_path, "_output_BFASTlite_breakpoint_SR_cal.gpkg", sep="_")
                               
                              st_write(breakpoint_SR, dsn = OutFile, layer = SRNames[layer])
                                
                              }

# Read the VALIDATION set of satellite data 
SR_GPKG = "./Intermediate product/cloud_free_product/__cloudfree_L8TS__SR_val.gpkg"
# Read each name layer of the satellite data
SRNames = st_layers(SR_GPKG)$name
SR = lapply(SRNames, function(name) st_read(SR_GPKG, layer=name))
# Remove geom column from the dataset
SR_nogeom <- lapply(SR, function(x) sf::st_drop_geometry(x))
# Run BFAST Lite on the dataset
for (layer in seq_along(SR)) { breakpoint_SR <- cal_BFAST("BIC", 0.25, response ~ trend, "trend", "trend",
                                                          SR_nogeom[[layer]], 
                                                          NULL,  plot = FALSE, cl = mycores)
                                # Write the output of BFAST lite
                                OutFile = paste(base_path, "_output_BFASTlite_breakpoint_SR_val.gpkg", sep="_")
                                
                                st_write(breakpoint_SR, dsn = OutFile, layer = SRNames[layer])

                              }
####################### METHOD 2 ###############################################
# Sum value of all band for each date. The output is one layer that contains summary values from all bands. Use this layer to run on the BFAST Lite.
# There is no need to calibrate in this method, but I run the calibration subdata set anyway
# Read the CALIBRATION set of satellite data 
SR_GPKG = "./Intermediate product/cloud_free_product/__cloudfree_L8TS__SR_cal.gpkg"
# Read each name layer of the satellite data
SRNames = st_layers(SR_GPKG)$name
SR = lapply(SRNames, function(name) st_read(SR_GPKG, layer=name))
# Convert a list of layers into a dataframe that contain all layers. Make sum of all rows for each sample id among all bands
SR_comb <- bind_rows(SR) %>%
  group_by( sample_id, centroid_x, centroid_y, tile, geom) %>%
  summarise_all(sum) %>%
  ungroup() %>%
  as.data.frame()


# Remove geom column from the dataset
SR_nogeom <- as.data.frame(SR_comb)
SR_nogeom <- SR_nogeom[, !names(SR_nogeom) %in% c( "geom")]
# Write the result
write.csv(SR_nogeom, file = "./Intermediate product/cloud_free_product/_cloudfree_L8TS_SR_comb_M2_cal.csv")
# Run the result on BFAST Lite
BFASTlite_BIC_T_025_SR_com_M2_cal <- cal_BFAST("BIC", 0.25, response ~ trend, "trend", "trend",
                                     "./Intermediate product/cloud_free_product/_cloudfree_L8TS_SR_comb_M2_cal.csv", 
                                     "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_T_025_SR_comb_M2_cal.csv",  plot = FALSE, cl = mycores)


# Read the VALIDATION dataset of satellite data
SR_GPKG = "./Intermediate product/cloud_free_product/__cloudfree_L8TS__SR_val.gpkg"
# Read each name layer of the satellite data
SRNames = st_layers(SR_GPKG)$name
SR = lapply(SRNames, function(name) st_read(SR_GPKG, layer=name))
# Convert a list of layers into a dataframe that contain all layers. Make sum of all rows for each sample id among all bands
SR_comb <- bind_rows(SR) %>%
  group_by( sample_id, centroid_x, centroid_y, tile, geom) %>%
  summarise_all(sum) %>%
  ungroup() %>%
  as.data.frame()

# Remove geom column from the dataset
SR_nogeom <- as.data.frame(SR_comb)
SR_nogeom <- SR_nogeom[, !names(SR_nogeom) %in% c("geom")]
# Write the result
write.csv(SR_nogeom, file = "./Intermediate product/cloud_free_product/_cloudfree_L8TS_SR_comb_M2_val.csv")
# Run the result on BFAST Lite
BFASTlite_BIC_T_025_SR_com_M2_val <- cal_BFAST("BIC", 0.25, response ~ trend, "trend", "trend",
                                               "./Intermediate product/cloud_free_product/_cloudfree_L8TS_SR_comb_M2_val.csv", 
                                               "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_T_025_SR_comb_M2_val.csv",  plot = FALSE, cl = mycores)



## TEST
# define paramaters for BFAST Monitor and preprocessing
params <- list(
  preprocessing = list(
    interpolate = FALSE,
    seasonality = ""),
  bfastLiteInputs = list(
    InputTS="",
    scrange=NULL, 
    scsig=0.05, 
    breaks="LWZ",
    sctype="OLS-MOSUM", 
    maginterval=0.1, 
    magcomponent="trend",
    magstat="RMSD", 
    magthreshold=-Inf, 
    coefcomponent="trend",
    coefthresholds=c(0, 0), 
    plot=FALSE, 
    quiet=TRUE, 
    order=3,
    formula=response ~ harmon, 
    TargetYears=NULL,
    seasonfreq=0.5, 
    breaknumthreshold=Inf, 
    altformula=NULL
  )
)

# load an preprocess the input datafile 
l8ndvi <- prepL8NDVIdataframe("./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg")

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
tslist <- pblapply(1:nrow(l8ndvi), 
                   FUN = getTrimmedts, 
                   dframe = l8ndvi, 
                   lookuptable = tab, 
                   tsfreq = params$preprocessing$seasonality, 
                   interpolate = params$preprocessing$interpolate, 
                   cl = mycores)

# --------- RUN BFAST Lite -------------------------------------------
runBFASTLite <- function(i, timeserieslist, parameters) {
  
  parameters$InputTS <- timeserieslist[[i]]
  return(do.call(UseBFASTLite, parameters))
}

print("Run BFAST Lite on every ts object")

bfastres <- lapply(1:length(tslist), runBFASTLite, 
                   timeserieslist = tslist, 
                   parameters = params$bfastLiteInputs)

NDVI_val <- read.csv("./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_T_025_val.csv")
