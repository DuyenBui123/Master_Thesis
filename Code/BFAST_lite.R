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
base_path <- "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\"
# load data ----------------------------------------------


cal_BFAST <- function( breaks, magthreshold, formula, magcomponent, coefcomponent, validating_data_path, output_path,plot, cl=mycores) {
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
      magcomponent=magcomponent,
      magstat="RMSD", 
      magthreshold= magthreshold, 
      coefcomponent=coefcomponent,
      coefthresholds=c(0, 0), 
      plot=plot, 
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
  if ("Magnitude.before" %in% colnames(theoutput)) {
    theoutput <- theoutput %>% 
      mutate(Magnitude = if_else(!is.na(Magnitude.RMSD),Magnitude.RMSD, 0))}
  
  # if (is.null(output_path) == FALSE) { # make an output file name 
  # outname <- output_path
  # 
  # 
  # # export
  # data.table::fwrite(theoutput,
  #                    file = outname,
  #                    showProgress = TRUE)
  # } 
  
  return(theoutput)
}


writeoutput <- function(outputbfast, ouputpath) {
  # join sample ids and centroid coordinates to the BFAST Monitor output
  print("Join output to sample ID information")
  theoutput <- pblapply(1:length(outputbfast), 
                        FUN = joinOutput, 
                        algo_list = outputbfast, 
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
  if ("Magnitude.before" %in% colnames(theoutput)) {
    theoutput <- theoutput %>% 
      mutate(Magnitude = if_else(!is.na(Magnitude.RMSD),Magnitude.RMSD, 0))}
  # make an output file name 
  outname <- ouputpath
  
  print(theoutput)
  # export
  data.table::fwrite(theoutput,
                     file = outname,
                     showProgress = TRUE)}
# Calibrate BIC, formula, order: 0-3, magnitude threshod: 0.25, 0.3, 1, -inf

# BIC, response~trend+harmon, order=0, magnitude threshold: 0.25



BFASTlite_BIC_SandT_025 <- cal_BFAST("BIC", 0.25, 
                                     response ~ harmon + trend, "trend", "trend", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_SandT_025.csv")



BFASTlite_BIC_T_025 <- cal_BFAST("BIC", 0.25, response ~ trend, "trend", "trend",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_T_025.csv")


BFASTlite_BIC_S_025 <- cal_BFAST("BIC", 0.25, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_S_025.csv")



BFASTlite_BIC_SandT_030 <- cal_BFAST("BIC", 0.30, 
                                     response ~ harmon + trend, "trend", "trend", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_SandT_030.csv")

BFASTlite_BIC_T_030 <- cal_BFAST("BIC", 0.30, response ~ trend, "trend", "trend",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_T_030.csv")


BFASTlite_BIC_S_030 <- cal_BFAST("BIC", 0.30, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_S_030.csv")

BFASTlite_BIC_SandT_050 <- cal_BFAST("BIC", 0.50, 
                                     response ~ harmon + trend, "trend", "trend", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_SandT_050.csv")

BFASTlite_BIC_T_050 <- cal_BFAST("BIC", 0.50, response ~ trend, "trend", "trend",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_T_050.csv")

BFASTlite_BIC_S_050 <- cal_BFAST("BIC", 0.50, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_S_050.csv")


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

# join sample ids and centroid coordinates to the BFAST Monitor output

BFASTlite_LWZ_SandT_025 <- cal_BFAST("LWZ", 0.25, 
                                     response ~ harmon + trend, "trend", "trend", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_SandT_025.csv", plot = FALSE, cl = mycores)

writeoutput(BFASTlite_LWZ_SandT_025, "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_SandT_025.csv")
BFASTlite_LWZ_T_025 <- cal_BFAST("LWZ", 0.25, response ~ trend, "trend", "trend",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_T_025.csv",plot = FALSE, cl = mycores)
writeoutput(BFASTlite_LWZ_T_025, "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_T_025.csv")
BFASTlite_LWZ_S_025 <- cal_BFAST("LWZ", 0.25, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_S_025.csv", plot = FALSE, cl = mycores)
writeoutput(BFASTlite_LWZ_S_025, "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_S_025.csv")
BFASTlite_LWZ_SandT_030 <- cal_BFAST("LWZ", 0.30, 
                                     response ~ harmon + trend, "trend", "trend", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_SandT_030.csv", plot = FALSE, cl = mycores)
writeoutput(BFASTlite_LWZ_SandT_030, "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_SandT_030.csv")
BFASTlite_LWZ_T_030 <- cal_BFAST("LWZ", 0.30, response ~ trend, "trend", "trend",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_T_030.csv", plot = FALSE, cl = mycores)
writeoutput(BFASTlite_LWZ_T_030, "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_T_030.csv")
BFASTlite_LWZ_S_030 <- cal_BFAST("LWZ", 0.30, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_S_030.csv", plot = FALSE, cl = mycores)
writeoutput(BFASTlite_LWZ_S_030, "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_S_030.csv")
BFASTlite_LWZ_SandT_050 <- cal_BFAST("LWZ", 0.50, 
                                     response ~ harmon + trend, "trend", "trend", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_SandT_050.csv", plot = FALSE, cl = mycores)
writeoutput(BFASTlite_LWZ_SandT_050, "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_SandT_050.csv")
BFASTlite_LWZ_T_050 <- cal_BFAST("LWZ", 0.50, response ~ trend, "trend", "trend",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_T_050.csv", plot = FALSE, cl = mycores)
writeoutput(BFASTlite_LWZ_T_050, "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_T_050.csv")
BFASTlite_LWZ_S_050 <- cal_BFAST("LWZ", 0.50, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_S_050.csv", plot = FALSE, cl = mycores)
writeoutput(BFASTlite_LWZ_S_050, "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_S_050.csv")
BFASTlite_LWZ_SandT_1 <- cal_BFAST("LWZ", 1, response ~ harmon + trend, "trend", "trend", 
                                   "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                   "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_SandT_1.csv", plot = FALSE, cl = mycores)
writeoutput(BFASTlite_LWZ_SandT_1, "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_SandT_1.csv")
BFASTlite_LWZ_T_1 <- cal_BFAST("LWZ", 1, response ~ trend, "trend", "trend",
                               "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                               "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_T_1.csv", plot = FALSE, cl = mycores)
writeoutput(BFASTlite_LWZ_T_1, "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_T_1.csv")
BFASTlite_LWZ_S_1 <- cal_BFAST("LWZ", 1, response ~ harmon,
                               c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                               c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                               "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg",
                               "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_S_1.csv", plot = FALSE, cl = mycores)
writeoutput(BFASTlite_LWZ_S_1, "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_S_1.csv")
BFASTlite_LWZ_SandT_inf <- cal_BFAST("LWZ", -Inf, response ~ harmon + trend, "trend", "trend", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_SandT_inf.csv", plot = FALSE, cl = mycores)
writeoutput(BFASTlite_LWZ_SandT_inf, "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_SandT_inf.csv")
BFASTlite_LWZ_T_inf <- cal_BFAST("LWZ", -Inf, response ~ trend, "trend", "trend",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg", 
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_T_inf.csv", plot = FALSE, cl = mycores)
writeoutput(BFASTlite_LWZ_T_inf, "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_T_inf.csv")
BFASTlite_LWZ_S_inf <- cal_BFAST("LWZ", -Inf, response ~ harmon,
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 c("harmoncos1","harmoncos2", "harmoncos3", "harmonsin1", "harmonsin2", "harmonsin3"),
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_cal.gpkg",
                                 "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_S_inf.csv", plot = FALSE, cl = mycores)
writeoutput(BFASTlite_LWZ_S_inf, "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_S_inf.csv")


###################################### WORK ON NDVI VALIDATION DATASET #################


BFASTlite_BIC_T_025_val <- cal_BFAST("BIC", 0.25, response ~ trend, "trend", "trend",
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_cloudfree_L8TS_NDVI_val.gpkg", 
                                     "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_T_025_val.csv",  plot = FALSE, cl = mycores)


###################################### RUN BFAST LITE ON SATELLITE DATA 

mainDir <- "C:\\Master_Thesis\\"
subDir <- "Output"
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)
SR_GPKG = "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\__cloudfree_L8TS__SR_val.gpkg"

SRNames = st_layers(SR_GPKG)$name
SR = lapply(SRNames, function(name) st_read(SR_GPKG, layer=name))
SR_nogeom <- lapply(SR, function(x) sf::st_drop_geometry(x))

for (layer in seq_along(SR)) { breakpoint_SR <- cal_BFAST("BIC", 0.25, response ~ trend, "trend", "trend",
                                                          SR_nogeom[[layer]], 
                                                          NULL,  plot = FALSE, cl = mycores)
OutFile = paste(base_path, "_output_BFASTlite_breakpoint_SR_.gpkg", sep="_")

st_write(breakpoint_SR, dsn = OutFile, layer = SRNames[layer])

}

