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
p_load(here, sf, zoo, data.table, tidyverse, pbapply,parallel, sandwich, hash) #, bfast, strucchangeRcpp)

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
mainDir <- "/home/duyen/Master_Thesis/"
subDir <- "Output"
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)
# load data ----------------------------------------------


# Calibrate BIC, formula: response ~ trend + harmon OR response ~ trend OR response ~ harmon, magnitude threshold: 0.25, 0.3, 1, -inf

BFASTlite_BIC_SandT_025 <- cal_BFAST("BIC", 0.25, 
                                     response ~ harmon + trend, "trend", "trend", 
                                     "./Intermediate product/cloud_free_product/BFAST_lite_NDVI_cal.gpkg", 
                                     "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_SandT_025.csv", plot = FALSE, cl=mycores)

##### adding date tiime
# theoutput <- data.table::rbindlist(BFASTlite_BIC_SandT_025$datetime, fill = TRUE)
# as.data.frame(BFASTlite_BIC_SandT_025$datetime[[1]])
# 
# datetime_dict <- hash()
# sampleid_NDVI_cal <- unique(BFASTlite_BIC_SandT_025$bp$sample_id)
# for (nr.sampid in 1: length(sampleid_NDVI_cal)) {
#   datetime_dict[sampleid_NDVI_cal[nr.sampid]] <- BFASTlite_BIC_SandT_025$datetime[[nr.sampid]]
# }

################################



BFASTlite_BIC_T_025 <- cal_BFAST("BIC", 0.25, response ~ trend, "trend", "trend",
                                 "./Intermediate product/cloud_free_product/BFAST_lite_NDVI_cal.gpkg", 
                                 "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_T_025.csv", plot = FALSE, cl=mycores)


BFASTlite_BIC_S_025 <- cal_BFAST("BIC", 0.25, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "./Intermediate product/cloud_free_product/BFAST_lite_NDVI_cal.gpkg",
                                 "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_S_025.csv", plot = FALSE, cl=mycores)


BFASTlite_BIC_SandT_050 <- cal_BFAST("BIC", 0.50, 
                                     response ~ harmon + trend, "trend", "trend", 
                                     "./Intermediate product/cloud_free_product/BFAST_lite_NDVI_cal.gpkg", 
                                     "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_SandT_050.csv", plot = FALSE, cl=mycores)

BFASTlite_BIC_T_050 <- cal_BFAST("BIC", 0.50, response ~ trend, "trend", "trend",
                                 "./Intermediate product/cloud_free_product/_cloudfree_L8TS_NDVI_cal.gpkg", 
                                 "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_T_050.csv", plot = FALSE, cl=mycores)

BFASTlite_BIC_S_050 <- cal_BFAST("BIC", 0.50, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "./Intermediate product/cloud_free_product/BFAST_lite_NDVI_cal.gpkg",
                                 "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_S_050.csv", plot = FALSE, cl=mycores)


BFASTlite_BIC_SandT_1 <- cal_BFAST("BIC", 1, response ~ harmon + trend, "trend", "trend",
                                   "./Intermediate product/cloud_free_product/BFAST_lite_NDVI_cal.gpkg", 
                                   "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_SandT_1.csv", plot = FALSE, cl=mycores)


BFASTlite_BIC_T_1 <- cal_BFAST("BIC", 1, response ~ trend, "trend", "trend",
                               "./Intermediate product/cloud_free_product/BFAST_lite_NDVI_cal.gpkg", 
                               "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_T_1.csv", plot = FALSE, cl = mycores)

BFASTlite_BIC_S_1 <- cal_BFAST("BIC", 1, response ~ harmon,
                               c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                               c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                               "./Intermediate product/cloud_free_product/BFAST_lite_NDVI_cal.gpkg",
                               "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_S_1.csv", plot = FALSE, cl = mycores)





# Calibrate LWZ, formula: response ~ trend + harmon OR response ~ trend OR response ~ harmon, magnitude threshold: 0.25, 0.3, 1, -inf

BFASTlite_LWZ_SandT_025 <- cal_BFAST("LWZ", 0.25, 
                                     response ~ harmon + trend, "trend", "trend", 
                                     "./Intermediate product/cloud_free_product/BFAST_lite_NDVI_cal.gpkg", 
                                     "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_SandT_025.csv", plot = FALSE, cl = mycores)


BFASTlite_LWZ_T_025 <- cal_BFAST("LWZ", 0.25, response ~ trend, "trend", "trend",
                                 "./Intermediate product/cloud_free_product/BFAST_lite_NDVI_cal.gpkg", 
                                 "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_T_025.csv",plot = FALSE, cl = mycores)

BFASTlite_LWZ_S_025 <- cal_BFAST("LWZ", 0.25, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "./Intermediate product/cloud_free_product/BFAST_lite_NDVI_cal.gpkg",
                                 "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_S_025.csv", plot = FALSE, cl = mycores)



BFASTlite_LWZ_SandT_050 <- cal_BFAST("LWZ", 0.50, 
                                     response ~ harmon + trend, "trend", "trend", 
                                     "./Intermediate product/cloud_free_product/BFAST_lite_NDVI_cal.gpkg", 
                                     "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_SandT_050.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_T_050 <- cal_BFAST("LWZ", 0.50, response ~ trend, "trend", "trend",
                                 "./Intermediate product/cloud_free_product/BFAST_lite_NDVI_cal.gpkg", 
                                 "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_T_050.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_S_050 <- cal_BFAST("LWZ", 0.50, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "./Intermediate product/cloud_free_product/BFAST_lite_NDVI_cal.gpkg",
                                 "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_S_050.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_SandT_1 <- cal_BFAST("LWZ", 1, response ~ harmon + trend, "trend", "trend", 
                                   "./Intermediate product/cloud_free_product/BFAST_lite_NDVI_cal.gpkg", 
                                   "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_SandT_1.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_T_1 <- cal_BFAST("LWZ", 1, response ~ trend, "trend", "trend",
                               "./Intermediate product/cloud_free_product/BFAST_lite_NDVI_cal.gpkg", 
                               "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_T_1.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_S_1 <- cal_BFAST("LWZ", 1, response ~ harmon,
                               c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                               c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                               "./Intermediate product/cloud_free_product/BFAST_lite_NDVI_cal.gpkg",
                               "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_S_1.csv", plot = FALSE, cl = mycores)
BFASTlite_LWZ_SandT_inf <- cal_BFAST("LWZ", -Inf, response ~ harmon + trend, "trend", "trend", 
                                     "./Intermediate product/cloud_free_product/BFAST_lite_NDVI_cal.gpkg", 
                                     "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_SandT_Inf.csv", plot = FALSE, cl = mycores)
BFASTlite_BIC_SandT_1 <- cal_BFAST("BIC", -Inf, response ~ harmon + trend, "trend", "trend", 
                                   "./Intermediate product/cloud_free_product/BFAST_lite_NDVI_cal.gpkg", 
                                   "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_SandT_Inf.csv", plot = FALSE, cl = mycores)


###################### Run BFAST lite ON NDVI VALIDATION DATASET ###############

# Run BFAST Lite on NDVI validation dataset
# Set of parameters is chosen based on F1 score of all calibration sets above. The one that has the highest F1 score is chosen
## RESULT: trend+harmon, LWZ, -Inf
BFASTlite_BIC_TandS_025_val <- cal_BFAST("LWZ", -Inf, response ~ trend + harmon, "trend", "trend",
                                         "./Intermediate product/cloud_free_product/BFAST_lite_NDVI_val.gpkg", 
                                         "./Intermediate product/cloud_free_product/_output_BFASTlite_LWZ_SandT_Inf_val.csv",  plot = FALSE, cl = mycores)

####################### METHOD 1 for running multiple bands of satellite data ##
# Select the most dominant floor year that contain more than or equal to 3 breakpoints 
# from all bands. Then choose the highest magnitude breakpoints.
###################### RUN BFAST LITE ON SATELLITE DATA  #######################
# Read the CALIBRATION set of satellite data 
SR_GPKG = "./Intermediate product/cloud_free_product/BFAST_lite_SR_cal.gpkg"
# Read each name layer of the satellite data
SRNames = st_layers(SR_GPKG)$name
SR = lapply(SRNames, function(name) st_read(SR_GPKG, layer=name))
# Remove geom column from the dataset
SR_nogeom <- lapply(SR, function(x) sf::st_drop_geometry(x))
# Initialize an empty list to store the lists from each iteration
all_lists <- list()
# Run BFAST Lite on the dataset

# for (layer in 2:7) { breakpoint_SR <- cal_BFAST("BIC", 0.25, response ~ trend + harmon, "trend", "trend",
#                                                           SR_nogeom[[layer]], 
#                                          NULL,  plot = FALSE, cl = mycores)
#                               # Write the output of BFAST lite
#                               OutFile = paste(base_path, "_output_BFASTlite_breakpoint_SR_cal.gpkg", sep="_")
#                                
#                               st_write(breakpoint_SR$bp, dsn = OutFile, layer = SRNames[layer])
#                               
#                              # Store the current list in all_lists
#                               all_lists[[paste0("list_", layer)]] <- breakpoint_SR$datetime
#                                 
# }

breakpoint_SR <- lapply(SR_nogeom, function(x) cal_BFAST("LWZ", -Inf, response ~ trend + harmon, "trend", "trend",
                                                         x, 
                                                         NULL,  plot = FALSE, cl = mycores))
for(layer in 2:7){
  # Write the output of BFAST lite
  OutFile = paste(base_path, "_output_BFASTlite_breakpoint_SR_cal.gpkg", sep="_")
  
  st_write(breakpoint_SR[[layer]]$bp, dsn = OutFile, layer = SRNames[layer])
  # Store the current list in all_lists
  all_lists[[paste0("list_", layer)]] <- breakpoint_SR[[layer]]$datetime
}


OutFile_datetime_SR_cal <- paste(base_path, "_output_BFASTlite_datetime_SR_cal.rds", sep="_")
# Save the overarching list into a file
saveRDS(all_lists, file = OutFile_datetime_SR_cal)


# Read the VALIDATION set of satellite data 
SR_GPKG = "./Intermediate product/cloud_free_product/BFAST_lite_SR_val.gpkg"
# Read each name layer of the satellite data
SRNames = st_layers(SR_GPKG)$name
SR = lapply(SRNames, function(name) st_read(SR_GPKG, layer=name))
# Remove geom column from the dataset
SR_nogeom <- lapply(SR, function(x) sf::st_drop_geometry(x))
all_lists <- list()
# Run BFAST Lite on the dataset
breakpoint_SR <- lapply(SR_nogeom[2:7], function(x) cal_BFAST("LWZ", -Inf, response ~ trend + harmon, "trend", "trend",
                                                              x, 
                                                              NULL,  plot = FALSE, cl = mycores))

for(layer in 1:6){
  # Write the output of BFAST lite
  OutFile = paste(base_path, "_output_BFASTlite_breakpoint_SR_val.gpkg", sep="_")
  
  st_write(breakpoint_SR[[layer]]$bp, dsn = OutFile, layer = SRNames[layer])
  # Store the current list in all_lists
  all_lists[[paste0("list_", layer)]] <- breakpoint_SR[[layer]]$datetime
}


OutFile_datetime_SR <- paste(base_path, "_output_BFASTlite_datetime_SR_val.rds", sep="_")
# Save the overarching list into a file
saveRDS(all_lists, file = OutFile_datetime_SR)
# ####################### METHOD 2 ###############################################
# # Sum value of all band for each date. The output is one layer that contains summary values from all bands. Use this layer to run on the BFAST Lite.
# # There is no need to calibrate in this method, but I run the calibration subdata set anyway
# # Read the CALIBRATION set of satellite data 
# SR_GPKG = "./Intermediate product/cloud_free_product/__cloudfree_L8TS__SR_cal.gpkg"
# # Read each name layer of the satellite data
# SRNames = st_layers(SR_GPKG)$name
# SR = lapply(SRNames, function(name) st_read(SR_GPKG, layer=name))
# # Convert a list of layers into a dataframe that contain all layers. Make sum of all rows for each sample id among all bands
# SR_comb <- bind_rows(SR) %>%
#   group_by( sample_id, centroid_x, centroid_y, tile, geom) %>%
#   summarise_all(sum) %>%
#   ungroup() %>%
#   as.data.frame()
# 
# 
# # Remove geom column from the dataset
# SR_nogeom <- as.data.frame(SR_comb)
# SR_nogeom <- SR_nogeom[, !names(SR_nogeom) %in% c( "geom")]
# # Write the result
# write.csv(SR_nogeom, file = "./Intermediate product/cloud_free_product/_cloudfree_L8TS_SR_comb_M2_cal.csv")
# # Run the result on BFAST Lite
# BFASTlite_BIC_T_025_SR_com_M2_cal <- cal_BFAST("BIC", 0.25, response ~ trend, "trend", "trend",
#                                      "./Intermediate product/cloud_free_product/_cloudfree_L8TS_SR_comb_M2_cal.csv", 
#                                      "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_T_025_SR_comb_M2_cal.csv",  plot = FALSE, cl = mycores)
# 
# 
# # Read the VALIDATION dataset of satellite data
# SR_GPKG = "./Intermediate product/cloud_free_product/__cloudfree_L8TS__SR_val.gpkg"
# # Read each name layer of the satellite data
# SRNames = st_layers(SR_GPKG)$name
# SR = lapply(SRNames, function(name) st_read(SR_GPKG, layer=name))
# # Convert a list of layers into a dataframe that contain all layers. Make sum of all rows for each sample id among all bands
# SR_comb <- bind_rows(SR) %>%
#   group_by( sample_id, centroid_x, centroid_y, tile, geom) %>%
#   summarise_all(sum) %>%
#   ungroup() %>%
#   as.data.frame()
# 
# # Remove geom column from the dataset
# SR_nogeom <- as.data.frame(SR_comb)
# SR_nogeom <- SR_nogeom[, !names(SR_nogeom) %in% c("geom")]
# # Write the result
# write.csv(SR_nogeom, file = "./Intermediate product/cloud_free_product/_cloudfree_L8TS_SR_comb_M2_val.csv")
# # Run the result on BFAST Lite
# BFASTlite_BIC_T_025_SR_com_M2_val <- cal_BFAST("BIC", 0.25, response ~ trend, "trend", "trend",
#                                                "./Intermediate product/cloud_free_product/_cloudfree_L8TS_SR_comb_M2_val.csv", 
#                                                "./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_T_025_SR_comb_M2_val.csv",  plot = FALSE, cl = mycores)
# 
# 

# ## TEST
# # define paramaters for BFAST Monitor and preprocessing
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
# 
# # load an preprocess the input datafile 
# l8ndvi_ <- prepL8NDVIdataframe("./Intermediate product/cloud_free_product/BFAST_lite_NDVI_val.gpkg")
# # load an preprocess the input datafile 
# l8ndvi_cal <- prepL8NDVIdataframe("./Intermediate product/cloud_free_product/BFAST_lite_NDVI_cal.gpkg")
# # l8ndvi_select <- l8ndvi[l8ndvi$sample_id %in% c(1042255, 1039666),]
# # inspect individual sample id 17873, 3495, 3510, 13990, 14043
# index_ <- which(l8ndvi_cal$sample_id ==1428507)
# # prepare time series ------------------------------------------------------
# 
# # get frequency table of acquisitions per year
# tab <- getFreqTabByYear("./Intermediate product/cloud_free_product/BFAST_lite_NDVI_val.gpkg")
# 
# # set the seasonality of the input time series
# if(!is.numeric(params$preprocessing$seasonality)){
#   # define the time series frequency by the minimum number of aquisitions in the time series
#   params$preprocessing$seasonality <- min(tab$Frequency[2:(length(tab$Frequency)-1)])
# }
# 
# # make list with timeseries (ts) objects 
# print("Prepare time series list")
# tslist_ <- pblapply(1:nrow(l8ndvi_),
#                    FUN = getTrimmedts, 
#                    dframe = l8ndvi_, 
#                    lookuptable = tab, 
#                    tsfreq = params$preprocessing$seasonality, 
#                    interpolate = params$preprocessing$interpolate, 
#                    cl = mycores)
# # inspect individual sample id
# tslist_se <- tslist_[index_:(index_+2)]
# 
# # --------- RUN BFAST Lite -------------------------------------------
# runBFASTLite <- function(i, timeserieslist, parameters) {
#   
#   parameters$InputTS <- timeserieslist[[i]]
#   return(do.call(UseBFASTLite, parameters))
# }
# 
# print("Run BFAST Lite on every ts object")
# 
# # bfastres_ <- lapply(1:500, runBFASTLite, 
# #                    timeserieslist = tslist, 
# #                    parameters = params$bfastLiteInputs)
# # inspect individual sample id
# bfastres_ <- lapply(1:2, runBFASTLite, 
#                                       timeserieslist = tslist_se, 
#                                        parameters = params$bfastLiteInputs)
# 
# theoutput <- pblapply(1:length(bfastres), 
#                       FUN = joinOutput, 
#                       algo_list = bfastres_, 
#                       dfids = l8ndvi[1:500,], 
#                       cl = mycores)
# theoutput <- data.table::rbindlist(theoutput, fill = TRUE)
# if ("Magnitude.before" %in% colnames(theoutput)) {
#   theoutput <- theoutput %>% 
#     mutate(Magnitude = if_else(!is.na(Magnitude.RMSD),Magnitude.RMSD, 0))}
# 
# outname <- "./Intermediate product/test_val.csv"
# 
# library(sf)
# # export
# data.table::fwrite(theoutput,
#                    file = outname,
#                    showProgress = TRUE)
# test_val <- read.csv("./Intermediate product/test_val.gpkg")
# test_val <- st_read("/home/duyen/Master_Thesis/Intermediate product/test_val.gpkg")
# NDVI_val <- read.csv("./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_SandT_025_val.csv")
# # Check file existence
# if (!file.exists("./Intermediate product/test_val.gpkg")) {
#   stop("File not found! Check the path.")
# }
# 
# # Read the GeoPackage
# test_val <- st_read("./Intermediate product/test_val.gpkg", quiet = FALSE)
# old_bfastres <- bfastres
# RMSD <- lapply(bfastres, function(x) return(x$Magnitude.RMSD[1]))
# 
# RMSD <- Filter(Negate(is.null), RMSD)
# 
# RMSE_cold <- lapply(bfastres, function(x) return(x$rmse[1]))
# RMSE_cold <- Filter(Negate(is.null), RMSE_cold)
# 
# RMSE_bfastlite <- lapply(bfastres, function(x) return(x$Magnitude.COLD_RMSE[1]))
# RMSE_bfastlite <- Filter(Negate(is.null), RMSE_bfastlite)
# x <- 1:64
# # Plot multiple lines using matplot
# matplot(x, cbind(RMSD, RMSE_cold, RMSE_bfastlite), type = "l", lty = 1, 
#         col = c("red", "blue", "green"), xlab = "X", 
#         ylab = "Y", main = "Multiple Lines Plot")
# legend("topright", legend = c("RMSD", "RMSE_cold", "RMSE_bfastlite"), 
#        col = c("red", "blue", "green"), 
#        lty = 1)
# runBFASTLite_datetime <- function(i, timeserieslist, parameters) {
#   
#   parameters$InputTS <- timeserieslist[[i]]
#   return(do.call(UseBFASTLite_datetime, parameters))
# }
# bfastres_datetime <- lapply(1:100, runBFASTLite_datetime, 
#                    timeserieslist = tslist, 
#                    parameters = params$bfastLiteInputs)
# 
# NDVI_val <- read.csv("./Intermediate product/cloud_free_product/_output_BFASTlite_BIC_T_025_val.csv")
