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
p_load(here, sf, zoo, data.table, tidyverse, pbapply,parallel, bfast, strucchangeRcpp)

# source external functions
here::i_am("Code\\BFAST_lite.R")
debugSource(here("GitHub_code", "Utils.R"))
source(here("GitHub_code", "getTrimmedts.R"))
source(here("GitHub_code", "BFASTutils.R"))
source(here("GitHub_code", "cglops-change-detection", "src", "bfast-cal", "04-validation.r"))

# add progress bar option to show it in the terminal 
pbo <- pboptions(type="timer")
mycores <- detectCores()
set.seed(123, kind = "L'Ecuyer-CMRG" )
# set pathes and create folder
base_path <- "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\"

# Create output folder
mainDir <- "C:\\Master_Thesis\\"
subDir <- "Output"
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)
# load data ----------------------------------------------


# Calibrate BIC, formula: response ~ trend + harmon OR response ~ trend OR response ~ harmon, magnitude threshold: 0.25, 0.3, 1, -inf

BFASTlite_BIC_SandT_025 <- cal_BFAST("BIC", 0.25, 
                                     response ~ harmon + trend, "trend", "trend", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_SandT_025.csv", plot = FALSE, cl=mycores)


BFASTlite_BIC_T_025 <- cal_BFAST("BIC", 0.25, response ~ trend, "trend", "trend",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_T_025.csv", plot = FALSE, cl=mycores)


BFASTlite_BIC_S_025 <- cal_BFAST("BIC", 0.25, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_S_025.csv", plot = FALSE, cl=mycores)



BFASTlite_BIC_SandT_030 <- cal_BFAST("BIC", 0.30, 
                                     response ~ harmon + trend, "trend", "trend", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_SandT_030.csv", plot = FALSE, cl=mycores)

BFASTlite_BIC_T_030 <- cal_BFAST("BIC", 0.30, response ~ trend, "trend", "trend",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_T_030.csv", plot = FALSE, cl=mycores)


BFASTlite_BIC_S_030 <- cal_BFAST("BIC", 0.30, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_S_030.csv", plot = FALSE, cl=mycores)

BFASTlite_BIC_SandT_050 <- cal_BFAST("BIC", 0.50, 
                                     response ~ harmon + trend, "trend", "trend", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_SandT_050.csv", plot = FALSE, cl=mycores)

BFASTlite_BIC_T_050 <- cal_BFAST("BIC", 0.50, response ~ trend, "trend", "trend",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_T_050.csv", plot = FALSE, cl=mycores)

BFASTlite_BIC_S_050 <- cal_BFAST("BIC", 0.50, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_S_050.csv", plot = FALSE, cl=mycores)


BFASTlite_BIC_SandT_1 <- cal_BFAST("BIC", 1, response ~ harmon + trend, "trend", "trend",
                                   "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                   "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_SandT_1.csv", plot = FALSE, cl=mycores)


BFASTlite_BIC_T_1 <- cal_BFAST("BIC", 1, response ~ trend, "trend", "trend",
                               "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                               "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_T_1.csv", plot = FALSE, cl = mycores)

BFASTlite_BIC_S_1 <- cal_BFAST("BIC", 1, response ~ harmon,
                               c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                               c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                               "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg",
                               "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_S_1.csv", plot = FALSE, cl = mycores)



BFASTlite_BIC_SandT_inf <- cal_BFAST("BIC", -Inf, response ~ harmon + trend, "trend", "trend",
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_SandT_Inf.csv", plot = FALSE, cl = mycores)

BFASTlite_BIC_T_inf <- cal_BFAST("BIC", -Inf, response ~ trend, "trend", "trend",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_T_Inf.csv", plot = FALSE, cl = mycores)


BFASTlite_BIC_S_inf <- cal_BFAST("BIC", -Inf, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_S_Inf.csv", plot = FALSE, cl = mycores)

# Calibrate LWZ, formula: response ~ trend + harmon OR response ~ trend OR response ~ harmon, magnitude threshold: 0.25, 0.3, 1, -inf

BFASTlite_LWZ_SandT_025 <- cal_BFAST("LWZ", 0.25, 
                                     response ~ harmon + trend, "trend", "trend", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_SandT_025.csv", plot = FALSE, cl = mycores)


BFASTlite_LWZ_T_025 <- cal_BFAST("LWZ", 0.25, response ~ trend, "trend", "trend",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_T_025.csv",plot = FALSE, cl = mycores)

BFASTlite_LWZ_S_025 <- cal_BFAST("LWZ", 0.25, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_S_025.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_SandT_030 <- cal_BFAST("LWZ", 0.30, 
                                     response ~ harmon + trend, "trend", "trend", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_SandT_030.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_T_030 <- cal_BFAST("LWZ", 0.30, response ~ trend, "trend", "trend",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_T_030.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_S_030 <- cal_BFAST("LWZ", 0.30, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_S_030.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_SandT_050 <- cal_BFAST("LWZ", 0.50, 
                                     response ~ harmon + trend, "trend", "trend", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_SandT_050.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_T_050 <- cal_BFAST("LWZ", 0.50, response ~ trend, "trend", "trend",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_T_050.csv", plot = TRUE, cl = mycores)

BFASTlite_LWZ_S_050 <- cal_BFAST("LWZ", 0.50, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_S_050.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_SandT_1 <- cal_BFAST("LWZ", 1, response ~ harmon + trend, "trend", "trend", 
                                   "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                   "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_SandT_1.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_T_1 <- cal_BFAST("LWZ", 1, response ~ trend, "trend", "trend",
                               "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                               "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_T_1.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_S_1 <- cal_BFAST("LWZ", 1, response ~ harmon,
                               c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                               c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                               "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg",
                               "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_S_1.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_SandT_inf <- cal_BFAST("LWZ", -Inf, response ~ harmon + trend, "trend", "trend", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_SandT_inf.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_T_inf <- cal_BFAST("LWZ", -Inf, response ~ trend, "trend", "trend",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_T_inf.csv", plot = FALSE, cl = mycores)

BFASTlite_LWZ_S_inf <- cal_BFAST("LWZ", -Inf, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_S_inf.csv", plot = FALSE, cl = mycores)



###################################### Run BFAST lite ON NDVI VALIDATION DATASET #################

# Run BFAST Lite on NDVI validation dataset
# Set of parameters is chosen based on F1 score of all calibration sets above. The one that has the highest F1 score is chosen
# BIC, response ~ trend, and magnitude threshold = 0.25 
BFASTlite_BIC_T_025_val <- cal_BFAST("BIC", 0.25, response ~ trend, "trend", "trend",
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_val.gpkg", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_T_025_val.csv",  plot = FALSE, cl = mycores)


###################################### RUN BFAST LITE ON SATELLITE DATA  #############################
# Read the CALIBRATION set of satellite data 
SR_GPKG = "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\__cloudfree_L8TS__SR_cal.gpkg"
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
"tile" %in% colnames(SR_nogeom[[1]])


# Read the VALIDATION set of satellite data 
SR_GPKG = "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\__cloudfree_L8TS__SR_val.gpkg"
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
############## METHOD 2 of aggregate the Satellite image
# Sum all bands
# Read the CALIBRATION set of satellite data 
SR_GPKG = "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\__cloudfree_L8TS__SR_cal.gpkg"
# Read each name layer of the satellite data
SRNames = st_layers(SR_GPKG)$name
SR = lapply(SRNames, function(name) st_read(SR_GPKG, layer=name))
# Remove geom column from the dataset
SR_nogeom <- lapply(SR, function(x) sf::st_drop_geometry(x))

