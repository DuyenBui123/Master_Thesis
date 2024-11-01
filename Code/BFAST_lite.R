install.packages("pacman")
library(pacman)
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(here, sf, zoo, data.table, tidyverse, pbapply,parallel, bfast, strucchangeRcpp)

# source external functions
here::i_am("Code\\BFAST_lite.R")
source(here("GitHub_code", "Utils.R"))
source(here("GitHub_code", "getTrimmedts.R"))
debugSource(here("GitHub_code", "BFASTutils.R"))
source(here("GitHub_code", "cglops-change-detection", "src", "bfast-cal", "04-validation.r"))

# add progress bar option to show it in the terminal 
pbo <- pboptions(type="timer")
mycores <- detectCores()
set.seed(123, kind = "L'Ecuyer-CMRG" )
datasetname <- "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\"
# load data ----------------------------------------------
debugSource(here("Code","BFAST_lite.R"))
l8ndvi <- prepL8NDVIdataframe("C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_L8_NDVI_output_BFASTLite.csv")
# define paramaters for BFAST Monitor and preprocessing

# CALIBRATE Breaks = LWZ, magthreshold = inf, formula = response ~ harmon + trend
params <- list(
  preprocessing = list(
    interpolate = FALSE,
    seasonality = ""),
  bfastLiteInputs = list(
    InputTS="",
    scrange=NULL, 
    scsig=0.05, 
    breaks= "LWZ",
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
    formula=response ~ harmon + trend, 
    TargetYears=NULL,
    seasonfreq=0.5, 
    breaknumthreshold=Inf, 
    altformula=NULL
  )
)
# load an preprocess the input datafile 
l8ndvi <- prepL8NDVIdataframe(here("Intermediate product", "cloud_free_product", "_cloudfree_L8TS_NDVI_cal.gpkg"))

# prepare time series ------------------------------------------------------

# get frequency table of acquisitions per year
tab <- getFreqTabByYear(here("Intermediate product", "cloud_free_product", "_cloudfree_L8TS_NDVI_cal.gpkg"))

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
runBFASTLite <- function(i,timeserieslist, parameters){
  parameters$InputTS <- timeserieslist[[i]]
  return(do.call(UseBFASTLite, parameters))
}

bfastres <- pblapply(1:length(tslist),
                       FUN = runBFASTLite,
                       timeserieslist = tslist,
                       parameters = params$bfastLiteInputs,
                       cl = mycores)

# remove now redundant tslist to clear memory
rm(tslist)
gc()

# join sample ids and centroid coordinates to the BFAST Monitor output
print("Join output to sample ID information")
theoutput <- pblapply(1:length(bfastres), 
                      FUN = joinOutput, 
                      algo_list = bfastres, 
                      dfids = l8ndvi, 
                      cl = mycores)
# concatenate all 1-row data.frames into one large data.frame
theoutput <- data.table::rbindlist(theoutput, fill = TRUE)

# add prediction ids to each prediction. 
theoutput <- theoutput %>% 
  group_by(sample_id) %>% 
  mutate(pred_id = row_number()*10) %>% 
  as.data.frame()



#' if there is a magnitude to report, 
#' fill the Magnitude column with "Mangitude.before"
theoutput <- theoutput %>% 
  mutate(Magnitude = if_else(!is.na(Magnitude.RMSD),Magnitude.RMSD, 0))
# make an output file name 
outname <- paste0(datasetname,'_L8_NDVI_output_BFASTLite.csv')


# export
data.table::fwrite(theoutput,
                   file = outname,
                   showProgress = TRUE)
# Validation
# CALIBRATE "LWZ", 0.25, response ~ harmon + trend


cal_BFAST <- function( breaks, magthreshold, formula, validating_data_path, output_path) {
params <- list(
  preprocessing = list(
    interpolate = FALSE,
    seasonality = ""),
  bfastLiteInputs = list(
    InputTS="",
    scrange=NULL, 
    scsig=0.05, 
    breaks= breaks,
    sctype="OLS-MOSUM", 
    maginterval=0.1, 
    magcomponent="trend",
    magstat="RMSD", 
    magthreshold= magthreshold, 
    coefcomponent="trend",
    coefthresholds=c(0, 0), 
    plot=FALSE, 
    quiet=TRUE, 
    order=3,
    formula=formula, 
    TargetYears=NULL,
    seasonfreq=0.5, 
    breaknumthreshold=Inf, 
    altformula=NULL
  )
)
# load an preprocess the input datafile 
l8ndvi <- prepL8NDVIdataframe(validating_data_path)

# prepare time series ------------------------------------------------------

# get frequency table of acquisitions per year
tab <- getFreqTabByYear(validating_data_path)

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
runBFASTLite <- function(i,timeserieslist, parameters){
  parameters$InputTS <- timeserieslist[[i]]
  return(do.call(UseBFASTLite, parameters))
}

bfastres <- pblapply(1:length(tslist),
                     FUN = runBFASTLite,
                     timeserieslist = tslist,
                     parameters = params$bfastLiteInputs,
                     cl = mycores)

# remove now redundant tslist to clear memory
rm(tslist)
gc()

# join sample ids and centroid coordinates to the BFAST Monitor output
print("Join output to sample ID information")
theoutput <- pblapply(1:length(bfastres), 
                      FUN = joinOutput, 
                      algo_list = bfastres, 
                      dfids = l8ndvi, 
                      cl = mycores)
# concatenate all 1-row data.frames into one large data.frame
theoutput <- data.table::rbindlist(theoutput, fill = TRUE)

# add prediction ids to each prediction. 
theoutput <- theoutput %>% 
  group_by(sample_id) %>% 
  mutate(pred_id = row_number()*10) %>% 
  as.data.frame()



#' if there is a magnitude to report, 
#' fill the Magnitude column with "Mangitude.before"
theoutput <- theoutput %>% 
  mutate(Magnitude = if_else(!is.na(Magnitude.RMSD),Magnitude.RMSD, 0))
# make an output file name 
outname <- output_path


# export
data.table::fwrite(theoutput,
                   file = outname,
                   showProgress = TRUE)
}

BFASTlite_BIC_SandT_025 <- cal_BFAST("BIC", 0.25, response ~ harmon + trend, "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_SandT_025.csv")
BFASTlite_BIC_T_025 <- cal_BFAST("BIC", 0.25, response ~ trend, "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_T_025.csv")
#######################################BFASTlite_BIC_S_025 <- cal_BFAST("BIC", 0.25, response ~ harmon, "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_S_025.csv")
params0 <- list(
  preprocessing = list(
    interpolate = FALSE,
    seasonality = ""),
  bfastLiteInputs = list(
    InputTS="",
    scrange=NULL, 
    scsig=0.05, 
    breaks= "BIC",
    sctype="OLS-MOSUM", 
    maginterval=0.1, 
    magcomponent=c("harmoncos1","harmoncos2"),
    magstat="RMSD", 
    magthreshold= 0.25, 
    coefcomponent=c("harmoncos1","harmoncos2"),
    coefthresholds=c(0, 0), 
    plot=FALSE, 
    quiet=TRUE, 
    order=3,
    formula= response ~ harmon, 
    TargetYears=NULL,
    seasonfreq=0.5, 
    breaknumthreshold=Inf, 
    altformula=NULL
  )
)
# load an preprocess the input datafile 
l8ndvi <- prepL8NDVIdataframe("C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg")

# prepare time series ------------------------------------------------------

# get frequency table of acquisitions per year
tab <- getFreqTabByYear("C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg")

# set the seasonality of the input time series

if(!is.numeric(params0$preprocessing$seasonality)){
  # define the time series frequency by the minimum number of aquisitions in the time series
  params0$preprocessing$seasonality <- min(tab$Frequency[2:(length(tab$Frequency)-1)])
}

# make list with timeseries (ts) objects 
print("Prepare time series list")
tslist <- pblapply(1:nrow(l8ndvi), 
                   FUN = getTrimmedts, 
                   dframe = l8ndvi, 
                   lookuptable = tab, 
                   tsfreq = params0$preprocessing$seasonality, 
                   interpolate = params0$preprocessing$interpolate, 
                   cl = mycores)
# --------- RUN BFAST Lite -------------------------------------------
runBFASTLite_BIC_025_S <- function(i,timeserieslist, parameters){
  parameters$InputTS <- timeserieslist[[i]]
  print(parameters)
  return(do.call(UseBFASTLite, parameters))
}

bfastres_BIC_025_S <- pblapply(1:length(tslist),
                     FUN = runBFASTLite_BIC_025_S,
                     timeserieslist = tslist,
                     parameters = params0$bfastLiteInputs,
                     cl = mycores)

# remove now redundant tslist to clear memory
rm(tslist)
gc()

# join sample ids and centroid coordinates to the BFAST Monitor output
print("Join output to sample ID information")
theoutput_BIC_025_S <- pblapply(1:length(bfastres_BIC_025_S), 
                              FUN = joinOutput, 
                              algo_list = bfastres_BIC_025_S, 
                              dfids = l8ndvi, 
                              cl = mycores)
# concatenate all 1-row data.frames into one large data.frame
theoutput_BIC_025_S <- data.table::rbindlist(theoutput_BIC_025_S, fill = TRUE)

# add prediction ids to each prediction. 
theoutput_BIC_025_S <- theoutput_BIC_025_S %>% 
  group_by(sample_id) %>% 
  mutate(pred_id = row_number()*10) %>% 
  as.data.frame()



#' if there is a magnitude to report, 
#' fill the Magnitude column with "Mangitude.before"
theoutput_BIC_025_S <- theoutput_BIC_025_S %>% 
  mutate(Magnitude = if_else(!is.na(Magnitude.RMSD),Magnitude.RMSD, 0))
# make an output file name 
outname <- "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_S_025.csv"


# export
data.table::fwrite(theoutput_BIC_025_S,
                   file = outname,
                   showProgress = TRUE)
BFASTlite_BIC_SandT_1 <- cal_BFAST("BIC", 1, response ~ harmon + trend, "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_SandT_1.csv")
#######################################BFASTlite_BIC_T_1 <- cal_BFAST("BIC", 1, response ~ trend, "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_T_1.csv")
params <- list(
  preprocessing = list(
    interpolate = FALSE,
    seasonality = ""),
  bfastLiteInputs = list(
    InputTS="",
    scrange=NULL, 
    scsig=0.05, 
    breaks= "BIC",
    sctype="OLS-MOSUM", 
    maginterval=0.1, 
    magcomponent="trend",
    magstat="RMSD", 
    magthreshold= 1, 
    coefcomponent="trend",
    coefthresholds=c(0, 0), 
    plot=FALSE, 
    quiet=TRUE, 
    order=3,
    formula=response ~ trend, 
    TargetYears=NULL,
    seasonfreq=0.5, 
    breaknumthreshold=Inf, 
    altformula=NULL
  )
)
# load an preprocess the input datafile 
l8ndvi <- prepL8NDVIdataframe("C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg")

# prepare time series ------------------------------------------------------

# get frequency table of acquisitions per year
tab <- getFreqTabByYear("C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg")

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
runBFASTLite_BIC_1_T <- function(i,timeserieslist, parameters){
  parameters$InputTS <- timeserieslist[[i]]
  return(do.call(UseBFASTLite, parameters))
}

bfastres_BIC_1_T <- pblapply(1:length(tslist),
                     FUN = runBFASTLite_BIC_1_T,
                     timeserieslist = tslist,
                     parameters = params$bfastLiteInputs,
                     cl = mycores)

# remove now redundant tslist to clear memory
rm(tslist)
gc()

# join sample ids and centroid coordinates to the BFAST Monitor output
print("Join output to sample ID information")
theoutput_BIC_1_T <- pblapply(1:length(bfastres_BIC_1_T), 
                              FUN = joinOutput, 
                              algo_list = bfastres_BIC_1_T, 
                              dfids = l8ndvi, 
                              cl = mycores)
# concatenate all 1-row data.frames into one large data.frame
theoutput_BIC_1_T <- data.table::rbindlist(theoutput_BIC_1_T, fill = TRUE)

# add prediction ids to each prediction. 
theoutput_BIC_1_T <- theoutput_BIC_1_T %>% 
  group_by(sample_id) %>% 
  mutate(pred_id = row_number()*10) %>% 
  as.data.frame()



#' if there is a magnitude to report, 
#' fill the Magnitude column with "Mangitude.before"
theoutput_BIC_1_T <- theoutput_BIC_1_T %>% 
  mutate(Magnitude = if_else(!is.na(Magnitude.RMSD),Magnitude.RMSD, 0))
# make an output file name 
outname <- "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_T_1.csv"


# export
data.table::fwrite(theoutput_BIC_1_T,
                   file = outname,
                   showProgress = TRUE)
BFASTlite_BIC_S_1 <- cal_BFAST("BIC", 1, response ~ harmon, "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_S_1.csv")
params <- list(
  preprocessing = list(
    interpolate = FALSE,
    seasonality = ""),
  bfastLiteInputs = list(
    InputTS="",
    scrange=NULL, 
    scsig=0.05, 
    breaks= "BIC",
    sctype="OLS-MOSUM", 
    maginterval=0.1, 
    magcomponent="trend",
    magstat="RMSD", 
    magthreshold= 1, 
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
l8ndvi <- prepL8NDVIdataframe("C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg")

# prepare time series ------------------------------------------------------

# get frequency table of acquisitions per year
tab <- getFreqTabByYear("C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg")

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
runBFASTLite_BIC_1_S <- function(i,timeserieslist, parameters){
  parameters$InputTS <- timeserieslist[[i]]
  return(do.call(UseBFASTLite, parameters))
}

bfastres_BIC_1_S <- pblapply(1:length(tslist),
                     FUN = runBFASTLite_BIC_1_S,
                     timeserieslist = tslist,
                     parameters = params$bfastLiteInputs,
                     cl = mycores)

# remove now redundant tslist to clear memory
rm(tslist)
gc()

# join sample ids and centroid coordinates to the BFAST Monitor output
print("Join output to sample ID information")
theoutput_BIC_1_S <- pblapply(1:length(bfastres_BIC_1_S), 
                      FUN = joinOutput, 
                      algo_list = bfastres_BIC_1_S, 
                      dfids = l8ndvi, 
                      cl = mycores)
# concatenate all 1-row data.frames into one large data.frame
theoutput_BIC_1_S <- data.table::rbindlist(theoutput_BIC_1_S, fill = TRUE)

# add prediction ids to each prediction. 
theoutput_BIC_1_S <- theoutput_BIC_1_S %>% 
  group_by(sample_id) %>% 
  mutate(pred_id = row_number()*10) %>% 
  as.data.frame()



#' if there is a magnitude to report, 
#' fill the Magnitude column with "Mangitude.before"
theoutput_BIC_1_S <- theoutput_BIC_1_S %>% 
  mutate(Magnitude = if_else(!is.na(Magnitude.RMSD),Magnitude.RMSD, 0))
# make an output file name 
outname <- "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_S_1.csv"


# export
data.table::fwrite(theoutput_BIC_1_S,
                   file = outname,
                   showProgress = TRUE)
BFASTlite_LWZ_SandT_1 <- cal_BFAST("LWZ", 1, response ~ harmon + trend, "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_SandT_1.csv")
BFASTlite_LWZ_T_1 <- cal_BFAST("LWZ", 1, response ~ trend, "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_T_1.csv")
params <- list(
  preprocessing = list(
    interpolate = FALSE,
    seasonality = ""),
  bfastLiteInputs = list(
    InputTS="",
    scrange=NULL, 
    scsig=0.05, 
    breaks= "LWZ",
    sctype="OLS-MOSUM", 
    maginterval=0.1, 
    magcomponent="trend",
    magstat="RMSD", 
    magthreshold= 1, 
    coefcomponent="trend",
    coefthresholds=c(0, 0), 
    plot=FALSE, 
    quiet=TRUE, 
    order=3,
    formula=response ~ trend, 
    TargetYears=NULL,
    seasonfreq=0.5, 
    breaknumthreshold=Inf, 
    altformula=NULL
  )
)
# load an preprocess the input datafile 
l8ndvi <- prepL8NDVIdataframe("C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg")

# prepare time series ------------------------------------------------------

# get frequency table of acquisitions per year
tab <- getFreqTabByYear("C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg")

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
runBFASTLite <- function(i,timeserieslist, parameters){
  parameters$InputTS <- timeserieslist[[i]]
  return(do.call(UseBFASTLite, parameters))
}

bfastres <- pblapply(1:length(tslist),
                     FUN = runBFASTLite,
                     timeserieslist = tslist,
                     parameters = params$bfastLiteInputs,
                     cl = mycores)

# remove now redundant tslist to clear memory
rm(tslist)
gc()

# join sample ids and centroid coordinates to the BFAST Monitor output
print("Join output to sample ID information")
theoutput <- pblapply(1:length(bfastres), 
                      FUN = joinOutput, 
                      algo_list = bfastres, 
                      dfids = l8ndvi, 
                      cl = mycores)
# concatenate all 1-row data.frames into one large data.frame
theoutput <- data.table::rbindlist(theoutput, fill = TRUE)

# add prediction ids to each prediction. 
theoutput <- theoutput %>% 
  group_by(sample_id) %>% 
  mutate(pred_id = row_number()*10) %>% 
  as.data.frame()



#' if there is a magnitude to report, 
#' fill the Magnitude column with "Mangitude.before"
theoutput <- theoutput %>% 
  mutate(Magnitude = if_else(!is.na(Magnitude.RMSD),Magnitude.RMSD, 0))
# make an output file name 
outname <- "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_T_1.csv"


# export
data.table::fwrite(theoutput,
                   file = outname,
                   showProgress = TRUE)
BFASTlite_LWZ_S_1 <- cal_BFAST("LWZ", 1, response ~ harmon, "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_S_1.csv")
params <- list(
  preprocessing = list(
    interpolate = FALSE,
    seasonality = ""),
  bfastLiteInputs = list(
    InputTS="",
    scrange=NULL, 
    scsig=0.05, 
    breaks= "LWZ",
    sctype="OLS-MOSUM", 
    maginterval=0.1, 
    magcomponent="trend",
    magstat="RMSD", 
    magthreshold= 1, 
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
l8ndvi <- prepL8NDVIdataframe("C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg")

# prepare time series ------------------------------------------------------

# get frequency table of acquisitions per year
tab <- getFreqTabByYear("C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg")

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
runBFASTLite <- function(i,timeserieslist, parameters){
  parameters$InputTS <- timeserieslist[[i]]
  return(do.call(UseBFASTLite, parameters))
}

bfastres <- pblapply(1:length(tslist),
                     FUN = runBFASTLite,
                     timeserieslist = tslist,
                     parameters = params$bfastLiteInputs,
                     cl = mycores)

# remove now redundant tslist to clear memory
rm(tslist)
gc()

# join sample ids and centroid coordinates to the BFAST Monitor output
print("Join output to sample ID information")
theoutput <- pblapply(1:length(bfastres), 
                      FUN = joinOutput, 
                      algo_list = bfastres, 
                      dfids = l8ndvi, 
                      cl = mycores)
# concatenate all 1-row data.frames into one large data.frame
theoutput <- data.table::rbindlist(theoutput, fill = TRUE)

# add prediction ids to each prediction. 
theoutput <- theoutput %>% 
  group_by(sample_id) %>% 
  mutate(pred_id = row_number()*10) %>% 
  as.data.frame()



#' if there is a magnitude to report, 
#' fill the Magnitude column with "Mangitude.before"
theoutput <- theoutput %>% 
  mutate(Magnitude = if_else(!is.na(Magnitude.RMSD),Magnitude.RMSD, 0))
# make an output file name 
outname <- "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_S_1.csv"


# export
data.table::fwrite(theoutput,
                   file = outname,
                   showProgress = TRUE)
