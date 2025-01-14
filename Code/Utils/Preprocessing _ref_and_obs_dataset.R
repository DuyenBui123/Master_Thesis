# ----  General info --------------------------
#' 
#' Name:             Duyen Bui
#' Study programs:   Master of Geo-information Science 
#' Student number:   1044767
#' 
#' Project:          Thesis: Compare the performance of multivariate and univariate time series models in detecting global land cover change
#' Supervisors:      Dainius Masiliunas & Nandin-Erdene Tsendbazar
#'
#' Script to preprocess reference data including define change or no change columns, merge WUR and IIASA data, calculate the statistics of change and no change based on observations and location id,
#' ploting land cover transition types.
#'
#'_____________________________________________________________________
#### Download and install packages ####
# LIBRARY
# create notin operator
`%notin%` <- Negate(`%in%`)

# Download packages if not available
pckgs <- c("tidyverse", "hash", "ggalluvial", "splitstackshape")

if (any(pckgs %notin% rownames(installed.packages())==TRUE)){
  install.packages(pckgs, repos = c(CRAN = "http://cloud.r-project.org"))}

# LOAD PACKAGE
sapply(pckgs, FUN = require, character.only = TRUE)
source("C:\\Master_Thesis\\Code\\Utils\\preprocessing_ref_and_sat_data.R")
library(ggalluvial)
library(hash)
library(splitstackshape)

# LOAD REFERECE DATA
WUR_data_path <- "C:\\Master_Thesis\\Data\\reference_global_100m_orig&change_year2015-2019_20210407.csv"
IIASA_data_path<-  "C:\\Master_Thesis\\Data\\Data_Global_quoted.csv"


# Read csv for reference data
WUR_data_raw <- read.csv(WUR_data_path, header=TRUE, stringsAsFactors=FALSE)
IIASA_data_raw <- read.csv(IIASA_data_path, header=TRUE, stringsAsFactors=FALSE)
IIASA_data <- read.csv("C:\\Master_Thesis\\Intermediate product\\IIASA_data.csv", header=TRUE, stringsAsFactors=FALSE) # This variable is the preprocessed data. Should use thid for analyzing to preventing preprocess the raw data again.
                                                                                                                # This data contains recalculated change, no change
                                                                                                                      # if run this line then dont run the two above
WUR_data <- read.csv("C:\\Master_Thesis\\Intermediate product\\WUR_data.csv", header=TRUE, stringsAsFactors=FALSE) # This variable is the preprocessed data. Should use thid for analyzing to preventing preprocess the raw data again.
# This data contains recalculated change, no change
# if run this line then dont run the two above

# PREPROCESSING REFERENCE DATA
# Select only necessary columns from the original dataset for WUR and IIASA dataset
WUR_data <- subset(WUR_data_raw, select = c(sample_id, subpix_mean_x, subpix_mean_y, location_id, bare, crops, grass, lichen, shrub, snow, trees, urban, water, fl.lichen, fl.grass, reference_year))
list_col_names_WUR <- c("sample_id", "centroid_x", "centroid_y", "location_id", "bare", "crops", "grass", "lichen", "shrub", "snow", "forest", "urban", "water", "lichen_wetland", "grass_wetland", "reference_year")
colnames(WUR_data) <- list_col_names_WUR

IIASA_data <- subset(IIASA_data_raw, select = c(sample_id,centroid_x, centroid_y, location_id, bare, crops, grassland, lichen_and_moss, shrub, snow_and_ice, tree, urban_built_up, water, wetland_herbaceous, reference_year))
list_col_names_IIASA <- c("sample_id", "centroid_x", "centroid_y", "location_id", "bare", "crops", "grass", "lichen", "shrub", "snow", "forest", "urban", "water", "wetland", "reference_year")
colnames(IIASA_data) <- list_col_names_IIASA
# Drop NA for IIASA
IIASA_data <- IIASA_data %>% drop_na()
summary(is.na(IIASA_data))
IIASA_data_raw <- IIASA_data_raw %>% drop_na()

# Redefine classes for WUR. Combine lichen- and grass- wetland into wetland
WUR_data$wetland <- WUR_data$lichen_wetland + WUR_data$grass_wetland
WUR_data <- subset(WUR_data, select = c(sample_id,centroid_x, centroid_y, location_id, bare, crops, grass, lichen, shrub, snow, forest, urban, water, wetland, reference_year))

# Create share variables for both dataset when preprocessing them
# Define land cover class included in two dataset
land_cover_type_WUR <- c("bare", "crops", "grass", "lichen", "shrub", "snow", "forest", "urban", "water", "wetland")
# Define the column numbers of the reference year and the first and the last land cover class of both dataset
refcol <- which(colnames(WUR_data)=="reference_year")
sample_id_col <- which(colnames(WUR_data)=="sample_id")
firstclass <- which(colnames(IIASA_data)=="bare")
lastclass <- which(colnames(IIASA_data)=="wetland")

############################### PREPROCESSING ####################################################################
#### WUR DATASET ####
# Check whether the fractional data need to be rescaled after adjusting land cover class columns
WUR_data$totalperc <- rowSums(WUR_data[,firstclass_WUR:lastclass_WUR])
if (length(unique(WUR_data$totalperc) != 1)) { print("The dataset need to be recaled")}
# Since the reference year of WUR dataset columns contain duplicated years, it needs to be adjusted to unique years
#Check reference year with only one year
one_ref_year <-0
for (sp_id in unique(WUR_data$sample_id)) {one_ref_year
  sp_id_wur <- WUR_data[WUR_data$sample_id == sp_id,]
  
  if(length(unique(sp_id_wur$reference_year))==1 & sp_id_wur[1,refcol] > 2015) {one_ref_year = one_ref_year + 1}}
one_ref_year

# If the reference year of a sample id contains only one year, change the reference years that are larger than 2015 to 2015
for (sp_id in unique(WUR_data$sample_id)) {
  sp_id_wur <- WUR_data[WUR_data$sample_id == sp_id,]
  
  if(length(unique(sp_id_wur$reference_year))==1 & sp_id_wur[1,refcol] > 2015) {WUR_data[WUR_data$sample_id == sp_id,]$reference_year <- 2015}}

# Check whether the sample ids that have only one reference year replaced to 2015
for (sp_id in unique(WUR_data$sample_id)) {one_ref_year
  sp_id_wur <- WUR_data[WUR_data$sample_id == sp_id,]
  
  if(length(unique(sp_id_wur$reference_year))==1 & sp_id_wur[1,refcol] > 2015) {print(sp_id_wur)}}

# OVERWRITE all duplicated years to the correct reference years

for (sp_id in unique(WUR_data$sample_id)) {sp_id_wur <- WUR_data[WUR_data$sample_id == sp_id,]
if(length(unique(sp_id_wur$reference_year))!=5) {sp_id_wur_dublicates <- sp_id_wur %>% 
  group_by(bare, crops, grass, lichen, shrub, snow, forest, urban, water, wetland, reference_year) %>%
  mutate(reference_year = case_when( n()>1 ~ reference_year + (row_number()-1),TRUE ~ reference_year)) %>%
  ungroup()
WUR_data[WUR_data$sample_id == sp_id,]$reference_year <- sp_id_wur_dublicates$reference_year}}

# Check whether all duplicated reference years are replaced

for (sp_id in unique(WUR_data$sample_id)) {sp_id_wur <- WUR_data[WUR_data$sample_id == sp_id,]
duplicated_ref_WUR <- duplicated(sp_id_wur[,sample_id_col:refcol])
if(any(duplicated_ref_WUR)==TRUE) {print(duplicated_ref_WUR)}}

for (sp_id in unique(WUR_data$sample_id)) {sp_id_wur <- WUR_data[WUR_data$sample_id == sp_id,]
duplicated_ref_WUR <- duplicated(sp_id_wur[,refcol])
if(any(duplicated_ref_WUR)==TRUE) {print(duplicated_ref_WUR)}}


#### Determine CHANGE OR NO CHANGE for each sample id ####

# Define year and number of reference year
year <- 2015
ref_nr_year <- 1
WUR_data$label <- c("no change")
# Make a for loop to run to all the unique sample id
for (sp_id in unique(WUR_data$sample_id)) {select_sp_id <- WUR_data[WUR_data$sample_id == sp_id,] 
year <- 2015
ref_nr_year <- 1
while (ref_nr_year < 5) {
  select_sp_id_1 <- select_sp_id[select_sp_id$reference_year == year,]
  select_sp_id_row1 <- select_sp_id_1 %>% select(bare, crops, grass, lichen, shrub, snow, forest, urban, water, wetland)
  
  select_sp_id_2 <- select_sp_id[select_sp_id$reference_year == year+1,]
  select_sp_id_row2 <- select_sp_id_2 %>% select(bare, crops, grass, lichen, shrub, snow, forest, urban, water, wetland)
  # Calculating change percentage of two consecutive rows
  if (length(select_sp_id_row1) == length(select_sp_id_row2)) {sub <- abs((select_sp_id_row1 - select_sp_id_row2))
  change_perc <- rowSums(sub, dims = 1)
  print(change_perc)
  ref_nr_year <- ref_nr_year+1
  # Adding threshold and define whether or not that year contain change, no change or small change
  if(change_perc >= 70) {WUR_data[WUR_data$sample_id == sp_id,][ref_nr_year,]$label <- "change"
  } else if (change_perc == 0) {WUR_data[WUR_data$sample_id == sp_id,][ref_nr_year,]$label <- "no change"
  } else {WUR_data[WUR_data$sample_id == sp_id,][ref_nr_year,]$label <- "small change"}
  year <- year+1

  } else{
  print("something wrong with sample id nr:", sp_id)
}}}

# WRITE a CSV file for WUR dataset contain change or no change labels and got preprocessed

write.csv(WUR_data,'C:\\Master_Thesis\\Intermediate product\\WUR_data.csv')
#### STATISTICS FOR CHANGE AND NO CHANGE for WUR ####
# Number of change
count_change <- 0
for (sp_id in unique(WUR_data$sample_id)) {
select_sp_id <- WUR_data[WUR_data$sample_id == sp_id,]
if((any(select_sp_id$label=="change") == TRUE)) {count_change <-count_change + 1
  
}}
count_change
# Number of small change
count_small_change <- 0
for (sp_id in unique(WUR_data$sample_id)) {
  select_sp_id <- WUR_data[WUR_data$sample_id == sp_id,]
  if((any(select_sp_id$label=="small change") == TRUE)) {count_small_change <-count_small_change + 1
  
  }}
count_small_change

# Number of no change
count_nochange <- 0
for (sp_id in unique(WUR_data$sample_id)) {
  select_sp_id <- WUR_data[WUR_data$sample_id == sp_id,]
  if((all(select_sp_id$label!="change") == TRUE)) {count_nochange <-count_nochange + 1
  
  }}
count_nochange
#_____________________________________________________________________________________________________________________________________________________________________#
#### IIASA DATASET ####
land_cover_type_IIASA <- c("bare", "crops", "grass", "lichen", "shrub", "snow", "forest", "urban", "water", "wetland")
firstclass_IIASA <- which(colnames(IIASA_data)=="bare")
lastclass_IIASA <- which(colnames(IIASA_data)=="wetland")

# Extract sample id that contain NaNs in the fractional class due to the removal of burnt, fallow crops, and not sure 
# and the values of these classes occupied 100 of the fractional values

remove_sampleid <-list() # Create an empty list for collecting sample id that contain 100% of burnt, or fallow_shifting_cultivation, or not_sure
for (sampleid in IIASA_data_raw$sample_id) {sample_id_IIASA <- IIASA_data_raw[IIASA_data_raw$sample_id == sampleid,]
if (any(sample_id_IIASA[["not_sure"]] >= 100)|any(sample_id_IIASA[["fallow_shifting_cultivation"]] >= 100)|any(sample_id_IIASA[["burnt"]] >= 100))
  # Check in the original dataset which sample id that have 100% of either "not sure", "fallow_shifting_cultivation", or "burnt"
{remove_sampleid <- append(remove_sampleid, sampleid)}} # Append the sample id in the remove list

remove_sampleid_IIASA <- unique(remove_sampleid) # Remove all the sample id contained in the list
extra_remove_sample_id <- c(1402822, 1404477) # There are two missing sample id that were not contained in the removal list 
#                                               due to I used drop NA values from IIASA_data_raw in the beginning, 
                                                  # thus some rows with 100% not sure or fallow are removed and left a few sample id
                                                  # with 3 or 4 rows leading to incomplete data
remove_sampleid_IIASA <- cbind(remove_sampleid_IIASA, extra_remove_sample_id)
IIASA_data <- IIASA_data[ ! IIASA_data$sample_id %in% remove_sampleid_IIASA, ] # A new IIASA dataset


# It is required to rescale the fractional data in IIASA.
# Rescaling
for (n in 1:nrow(IIASA_data)) {
  totalvalue <- rowSums(IIASA_data[n,][,firstclass_IIASA:lastclass_IIASA])
  if(isTRUE(totalvalue != 100)==TRUE){
    for (classtype in land_cover_type_IIASA) {
      IIASA_data[n,][[classtype]] <- (IIASA_data[n,][[classtype]]/totalvalue)*100
    }
  }}
# Check for NANs
IIASA_data$totalperc <- rowSums(IIASA_data[,firstclass_IIASA:lastclass_IIASA])
for (n in 1:nrow(IIASA_data)) {
  
  if(is.nan(IIASA_data$totalperc[n]) == TRUE){print(IIASA_data[n,])}}
# Check if there is any duplicated reference year

count_dublicated_row <- 0
for (sp_id in unique(IIASA_data$sample_id)) {sp_id_IIASA <- IIASA_data[IIASA_data$sample_id == sp_id,]
duplicated_ref_IIASA <- duplicated(sp_id_IIASA[,sample_id_col:refcol])
if(any(duplicated_ref_IIASA)==TRUE) {count_dublicated_row <- count_dublicated_row + 1
                                            print (sp_id)
                                            print(duplicated_ref_IIASA)}}


####################### TEST for change or no change function ##################
change_nochange_label(IIASA_data,IIASA_data,folder = "C:/Master_Thesis/Intermediate product/")


#############################################################################
#### Determine CHANGE OR NO CHANGE for each sample id ####

## Define year and number of reference year
year <- 2015
ref_nr_year <- 1
IIASA_data$label <- c("no change")

## Make a for loop to run to all the unique sample id
options(warn=1)
for (sp_id in unique(IIASA_data$sample_id)) {select_sp_id <- IIASA_data[IIASA_data$sample_id == sp_id,] 
year <- 2015
ref_nr_year <- 1
while (ref_nr_year < 4) {
  select_sp_id_1 <- select_sp_id[select_sp_id$reference_year == year,]
  select_sp_id_row1 <- select_sp_id_1 %>% select(bare, crops, grass, lichen, shrub, snow, forest, urban, water, wetland)
  
  select_sp_id_2 <- select_sp_id[select_sp_id$reference_year == year+1,]
  select_sp_id_row2 <- select_sp_id_2 %>% select(bare, crops, grass, lichen, shrub, snow, forest, urban, water, wetland)
  # Calculating change percentage for two consecutive rows
  sub <- abs((select_sp_id_row1 - select_sp_id_row2))
  change_perc <- rowSums(sub, dims = 1)
  print(change_perc)
  ref_nr_year <- ref_nr_year+1
  # Adding threshold to determine change labels
  if(change_perc >= 70) {IIASA_data[IIASA_data$sample_id == sp_id,][ref_nr_year,]$label <- "change"
  } else if (change_perc == 0) {IIASA_data[IIASA_data$sample_id == sp_id,][ref_nr_year,]$label <- "no change"
  } else {IIASA_data[IIASA_data$sample_id == sp_id,][ref_nr_year,]$label <- "small change"}
  year <- year+1
  print(sp_id)
}}
warnings()
# Write a csv for preprocessed and labeled IIASA dataset
write.csv(IIASA_data,'C:/Master_Thesis/Intermediate product/IIASA_data.csv', row.names = FALSE)

#### STATISTICS FOR CHANGE AND NO CHANGE for IIASA ####
# Number of change
count_change <- 0
for (sp_id in unique(IIASA_data$sample_id)) {
  select_sp_id <- IIASA_data[IIASA_data$sample_id == sp_id,]
  if((any(select_sp_id$label=="change") == TRUE)) {count_change <-count_change + 1
  
  }}
count_change
# Number of small change
count_small_change <- 0
for (sp_id in unique(IIASA_data$sample_id)) {
  select_sp_id <- IIASA_data[IIASA_data$sample_id == sp_id,]
  if((any(select_sp_id$label=="small change") == TRUE)) {count_small_change <-count_small_change + 1
  
  }}
count_small_change
# Number of no change
count_nochange <- 0
for (sp_id in unique(IIASA_data$sample_id)) {
  select_sp_id <- IIASA_data[IIASA_data$sample_id == sp_id,]
  if((all(select_sp_id$label!="change") == TRUE)) {count_nochange <-count_nochange + 1
  
  }}
count_nochange

#__________________________________________________________________________________________________________________________________________________#
#### COMBINE WUR AND IIASA ####
# Make IIASA and WUR data matched
# Select columns and rearrange the columns in IIASA to make it match with the format and order in WUR data
IIASA_data <- IIASA_data %>% select(X, sample_id , centroid_x ,centroid_y ,location_id , bare ,crops ,grass ,lichen ,shrub , snow , forest ,urban , water ,wetland , reference_year     , label ,reclass)
# Combine two reference dataset
com_ref_data <- rbind(WUR_data, IIASA_data)
# Create a dictionary to describe which sample id contain change and which not.
# This will be helpful in spliting calibration and validation data based on sample id 
# and make sure that the amount of change and no change is distributed equally
unique_sampleid <- list()
for (sampleid in unique(com_ref_data$sample_id)) {sampleid <- as.character(sampleid) # Create a list of unique sample id
  unique_sampleid <- append(unique_sampleid, sampleid)}

unique_sampleid_dict <- hash() # Create a dictionary that contains unique sample id as keys and change or no change as value

for (sampleid in unique_sampleid) {if (!has.key(sampleid, unique_sampleid_dict)) {#https://stackoverflow.com/questions/39492564/how-to-create-a-dictionary-and-insert-key-with-a-list-of-value-in-r
  unique_sampleid_dict[sampleid] <- "no status"}} # First set no status to all rows

for (sampleid in unique(com_ref_data$sample_id)){ sample_id <- com_ref_data[com_ref_data$sample_id ==  sampleid,]
  if (any(sample_id$label == "change")) { sampleid <- as.character(sampleid)
  unique_sampleid_dict[[sampleid]] <- "change"} else {sampleid <- as.character(sampleid)
                                                       unique_sampleid_dict[[sampleid]] <- "no change"}
}

# Convert the dictionary to a data frame
# Extract keys and values from the hash object
keys <- keys(unique_sampleid_dict)
values <- values(unique_sampleid_dict)

# Create a data frame from the keys and values
unique_sampleid_df <- data.frame(
  sample_id = keys,
  label = values, 
  check.rows = TRUE
)

#### SPLIT CALIBRATE DATA AND VALIDATA DATA ####
set.seed(1603)
# Using stratified random method to split the data. Calibrated data has 30% of the dataset and validated data contains 70%
cal_data <- stratified(unique_sampleid_df, c("label"), .3, select = list(label = c("change", "no change")))
val_data <- subset(unique_sampleid_df, !(unique_sampleid_df$sample_id %in% cal_data$sample_id))
cal_data_int <- as.integer(cal_data$sample_id)
val_data_int <- as.integer(val_data$sample_id)
# Apply the split from the dictionary to the real data
# Split calibration reference data
cal_compl_data <- subset(com_ref_data, (com_ref_data$sample_id %in% cal_data_int)) 
# Add "change_at_100m" and "change_yr" columns
cal_compl_data <- cal_compl_data %>% mutate(change_yr = as.double(reference_year))

for (i in 1:length(cal_compl_data$change_at_100m)) {
  if (cal_compl_data$change_at_100m[i] == "FALSE") {
    cal_compl_data$change_yr[i] <- ""
  } else {
    cal_compl_data$change_yr[i] <- cal_compl_data$reference_year[i]
  }
}

# Convert integer to numeric with three decimal for change_yr column
cal_compl_data$change_yr <- noquote(format(cal_compl_data$change_yr, digits=3, nsmall=3))
cal_compl_data$change_at_100m <- cal_compl_data$label
# Change the value in change_at_100m from change to TRUE, no change and small change to FALSE
cal_compl_data$change_at_100m <- replace(cal_compl_data$change_at_100m , cal_compl_data$change_at_100m == "change"  , "TRUE")
cal_compl_data$change_at_100m <- replace(cal_compl_data$change_at_100m , cal_compl_data$change_at_100m == "no change"|cal_compl_data$change_at_100m == "small change", "FALSE")
# Write a csv file for calibration data
file_path_cal <- "C:/Master_Thesis/Intermediate product/cal_compl_data.csv"
write.csv(cal_compl_data, file_path_cal)
# Split validation reference data
val_compl_data <- subset(com_ref_data, (com_ref_data$sample_id %in% val_data_int)) # Create validation data
# Add "change_at_100m" and "change_yr" columns
val_compl_data <- val_compl_data %>% mutate(change_yr = as.double(reference_year))

for (i in 1:length(val_compl_data$change_at_100m)) {
  if (val_compl_data$change_at_100m[i] == "FALSE") {
    val_compl_data$change_yr[i] <- ""
  } else {
    val_compl_data$change_yr[i] <- val_compl_data$reference_year[i]
  }
}
# Convert integer to numeric with three decimal for change_yr column
val_compl_data$change_yr <- noquote(format(val_compl_data$change_yr, digits=3, nsmall=3))
val_compl_data$change_at_100m <- val_compl_data$label
# Change the value in change_at_100m from change to TRUE, no change and small change to FALSE
val_compl_data$change_at_100m <- replace(val_compl_data$change_at_100m , val_compl_data$change_at_100m == "change"  , "TRUE")
val_compl_data$change_at_100m <- replace(val_compl_data$change_at_100m , val_compl_data$change_at_100m == "no change"|val_compl_data$change_at_100m == "small change", "FALSE")

# Write a csv file for validation data
write.csv(val_compl_data,'C:\\Master_Thesis\\Intermediate product\\val_compl_data.csv')

#### STATISTICS FOR CHANGE AND NO CHANGE for combined data ####

change_nochange_based_sampleid(cal_compl_data)

change_nochange_based_sampleid(val_compl_data)

# Calculate number of unique locations that contain each land cover class
land_cover_data <- list()
for(name_class in land_cover_type_WUR){
  land_cover_data[[paste0(name_class, "_data")]] <- list()
}

for (sp_id in unique(com_ref_data$sample_id)) {
  select_sp_id <- com_ref_data[com_ref_data$sample_id == sp_id,]
  for (name_class in land_cover_type_WUR) { 
    if(any(select_sp_id[[name_class]] != 0))
    {land_cover_data[[paste0(name_class, "_data")]] <- append(land_cover_data[[paste0(name_class, "_data")]], sp_id)
    }
  }
}
# Print out the length of each class
for(name_class in land_cover_type_WUR){
  variablename <- paste0(name_class, "_data")
  print(paste("length of", variablename,":", length(land_cover_data[[variablename]])))
}
# Calculate the number of unique locations that contains change in each land cover class for calibration data
land_cover_data_com <- list()
for(name_class in land_cover_type_WUR){
  land_cover_data_com[[paste0(name_class, "_data")]] <- list()
}

for (sp_id in unique(com_ref_data$sample_id)) {
  select_sp_id <- com_ref_data[com_ref_data$sample_id == sp_id,]
  for (name_class in land_cover_type_WUR) { 
    if(any(select_sp_id[[name_class]] != 0) & (any(select_sp_id$label == "change")))
    {land_cover_data_com[[paste0(name_class, "_data")]] <- append(land_cover_data_com[[paste0(name_class, "_data")]], sp_id)
    }
  }
}
# Print out the length of each land cover class that contains change
for(name_class in land_cover_type_WUR){
  variablename <- paste0(name_class, "_data")
  print(paste("length of", variablename,":", length(land_cover_data_com[[variablename]])))
}
# Calculate the number of unique locations that contains change in each land cover class for validation data
land_cover_data_val <- list()
for(name_class in land_cover_type_WUR){
  land_cover_data_val[[paste0(name_class, "_data")]] <- list()
}

for (sp_id in unique(val_compl_data$sample_id)) {
  select_sp_id <- val_compl_data[val_compl_data$sample_id == sp_id,]
  for (name_class in land_cover_type_WUR) { 
    if(any(select_sp_id[[name_class]] != 0) & (any(select_sp_id$label == "change")))
    {land_cover_data_val[[paste0(name_class, "_data")]] <- append(land_cover_data_val[[paste0(name_class, "_data")]], sp_id)
    }
  }
}
for(name_class in land_cover_type_WUR){
  variablename <- paste0(name_class, "_data")
  print(paste("length of", variablename,":", length(land_cover_data_val[[variablename]])))}

#### MAKING ALLUVIA DIAGRAM ####
# Creating a dictionaries to store types of land cover transition types and the frequency of transition happened in the corresponding transition for calibration data
cal_transition_types <- from_to_transition_type(cal_compl_data)
# Write csvs for each group of output: 1, all transitions including non transition (all);
# 2, only transitions (only); 3, only transition types that are interested in this study selected (interested)
write.csv(cal_transition_types$all,'C:\\Master_Thesis\\Intermediate product\\cal_all_transition_types.csv')
write.csv(cal_transition_types$only,'C:\\Master_Thesis\\Intermediate product\\cal_only_transition_types.csv')
write.csv(cal_transition_types$interested,'C:\\Master_Thesis\\Intermediate product\\cal_interested_transition_types.csv')
# Creating a dictionaries to store types of land cover transition types and the frequency of transition happened in the corresponding transition for validation data
val_transition_types <- from_to_transition_type(val_compl_data)
# Write csvs for each group of output: 1, all transitions including non transition (all);
# 2, only transitions (only); 3, only transition types that are interested in this study selected (interested)
write.csv(val_transition_types$all,'C:\\Master_Thesis\\Intermediate product\\val_all_transition_types.csv')
write.csv(val_transition_types$only,'C:\\Master_Thesis\\Intermediate product\\val_only_transition_types.csv')
write.csv(val_transition_types$interested,'C:\\Master_Thesis\\Intermediate product\\val_interested_transition_types.csv')

#############################################################################################################################################################
# STatistics for transition type
cal_all_transition_types <- read.csv('C:\\Master_Thesis\\Intermediate product\\cal_all_transition_types.csv', header = TRUE)
# Calculate the total number of a land cover class transit to other classes
cal_all_transition_types <- cal_all_transition_types %>% mutate(from_count = sum(Count), .by = From)
# Calculate the total number of each land cover class that was converted into from the "From" columns
cal_all_transition_types <- cal_all_transition_types %>% mutate(to_count = sum(Count), .by = To)

cal_only_transition_types <- read.csv('C:\\Master_Thesis\\Intermediate product\\cal_only_transition_types.csv')
# Calculate the total number of a land cover class transit to other classes
cal_only_transition_types <- cal_only_transition_types %>% mutate(from_count = sum(Count), .by = From)
# Calculate the total number of each land cover class that was converted into from the "From" columns
cal_only_transition_types <- cal_only_transition_types %>% mutate(to_count = sum(Count), .by = To)

cal_interested_transition_types <- read.csv('C:\\Master_Thesis\\Intermediate product\\cal_interested_transition_types.csv')
# Calculate the total number of a land cover class transit to other classes
cal_interested_transition_types <- cal_interested_transition_types %>% mutate(from_count = sum(Count), .by = From)
# Calculate the total number of each land cover class that was converted into from the "From" columns
cal_interested_transition_types <- cal_interested_transition_types %>% mutate(to_count = sum(Count), .by = To)


val_all_transition_types <- read.csv('C:\\Master_Thesis\\Intermediate product\\val_all_transition_types.csv', header = TRUE)
# Calculate the total number of a land cover class transit to other classes
val_all_transition_types <- val_all_transition_types %>% mutate(from_count = sum(Count), .by = From)
# Calculate the total number of each land cover class that was converted into from the "From" columns
val_all_transition_types <- val_all_transition_types %>% mutate(to_count = sum(Count), .by = To)

val_only_transition_types <- read.csv('C:\\Master_Thesis\\Intermediate product\\val_only_transition_types.csv')
# Calculate the total number of a land cover class transit to other classes
val_only_transition_types <- val_only_transition_types %>% mutate(from_count = sum(Count), .by = From)
# Calculate the total number of each land cover class that was converted into from the "From" columns
val_only_transition_types <- val_only_transition_types %>% mutate(to_count = sum(Count), .by = To)

val_interested_transition_types <- read.csv('C:\\Master_Thesis\\Intermediate product\\val_interested_transition_types.csv')
# Calculate the total number of a land cover class transit to other classes
val_interested_transition_types <- val_interested_transition_types %>% mutate(from_count = sum(Count), .by = From)
# Calculate the total number of each land cover class that was converted into from the "From" columns
val_interested_transition_types <- val_interested_transition_types %>% mutate(to_count = sum(Count), .by = To)

# Making alluvia diagram
# For calibration

# For all
# Specify colors for each class
cols <- c("red", "yellow", "darkgreen", "green", "deeppink", "lightblue", "#FFA500", "#00FFFF", "blue", "#698B22")

ggplot(as.data.frame(cal_all_transition_types),
       aes(y = Count, axis1 = From, axis2 = To)) +
  geom_alluvium(alpha = 1, aes(fill = From), width = 1/12) # Specify the transparency, the width and the color of alluvium
+ geom_stratum(stat = "stratum", width = 1/10, aes(fill = after_stat(stratum))) +
  #geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  geom_text(stat = "stratum", aes(label= cal_all_transition_types$from_count, position = "dodge")) # Specify the position of the text
+ geom_text(stat = "stratum", aes(label= cal_all_transition_types$to_count)) +
  scale_x_discrete(limits = c("From", "To"), expand = c(.05, .05)) +
  
  scale_color_manual(values=cols, aesthetics = c("colour", "fill")) +
  ggtitle("All land cover transition - calibration data") + theme(plot.title = element_text(size=22)) + theme(plot.title=element_text(hjust=0.5))
  theme_void()
  ggsave("All land cover transition - calibration data.pdf", path = "C:\\Master_Thesis\\Output", scale = 2)
# For only
  cols <- c("red", "yellow", "darkgreen", "green", "deeppink", "lightblue", "#FFA500", "#00FFFF", "blue", "#698B22")
  
  ggplot(as.data.frame(cal_only_transition_types),
         aes(y = Count, axis1 = From, axis2 = To)) +
    geom_alluvium(alpha = 1, aes(fill = From), width = 1/12) +
    
    geom_stratum(stat = "stratum", width = 1/10, aes(fill = after_stat(stratum))) +
    #geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    geom_text(stat = "stratum", aes(label= cal_only_transition_types$from_count, position = "dodge")) +
    geom_text(stat = "stratum", aes(label= cal_only_transition_types$to_count)) +
    scale_x_discrete(limits = c("From", "To"), expand = c(.05, .05)) +
    
    scale_color_manual(values=cols, aesthetics = c("colour", "fill")) +
    ggtitle("Land cover transition only - calibration data") + theme(plot.title = element_text(size=22)) + theme(plot.title=element_text(hjust=0.5))
  theme_void()
  ggsave("Land cover transition only_calibration data.pdf", path = "C:\\Master_Thesis\\Output", scale = 2)

# For interested 
  cols <- c("red", "darkgreen", "green", "#00FFFF", "blue", "#698B22", "yellow", "lightblue", "#FFA500")
  
  ggplot(as.data.frame(cal_interested_transition_types),
         aes(y = Count, axis1 = From, axis2 = To)) +
    geom_alluvium(alpha = 1, aes(fill = From), width =  1/10, curve_type = "arctangent") +
    
    geom_stratum(stat = "stratum", width = 1/10, aes(fill = after_stat(stratum))) +
    #geom_label(stat = "stratum", aes(label = after_stat(stratum)), how.legend = FALSE, vjust = -0.5) +
    geom_text(stat = "stratum", aes(label= cal_interested_transition_types$from_count, position = "dodge")) +
    geom_text(stat = "stratum", aes(label= cal_interested_transition_types$to_count)) + #https://stackoverflow.com/questions/60725388/how-to-add-percentage-values-to-strata-in-alluvial-plot-with-ggalluvial
    scale_x_discrete(limits = c("From", "To"), expand = c(.05, .05)) +
    scale_color_manual(values=cols, aesthetics = c("colour", "fill")) +
    ggtitle("Interested land cover transition only - calibration data") + theme(plot.title = element_text(size=22)) + theme(plot.title=element_text(hjust=0.5))
  theme_void()
  ggsave("Interested land cover transition only - calibration data.pdf", path = "C:\\Master_Thesis\\Output", scale = 2)
  
# FOR VALIDATAION

  cols <- c("red", "yellow", "darkgreen", "green", "deeppink", "lightblue", "#FFA500", "#00FFFF", "blue", "#698B22")
  
  ggplot(as.data.frame(val_all_transition_types),
         aes(y = Count, axis1 = From, axis2 = To)) +
    geom_alluvium(alpha = 1, aes(fill = From), width = 1/12) +
    
    geom_stratum(stat = "stratum", width = 1/10, aes(fill = after_stat(stratum))) +
    #geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    geom_text(stat = "stratum", aes(label= val_all_transition_types$from_count, position = "dodge")) +
    geom_text(stat = "stratum", aes(label= val_all_transition_types$to_count)) +
    scale_x_discrete(limits = c("From", "To"), expand = c(.05, .05)) +
    
    scale_color_manual(values=cols, aesthetics = c("colour", "fill")) +
    ggtitle("All land cover transition - validation data") + theme(plot.title = element_text(size=22)) + theme(plot.title=element_text(hjust=0.5))
  theme_void()
  ggsave("All land cover transition - validation data.pdf", path = "C:\\Master_Thesis\\Output", scale = 2)


  cols <- c("red", "yellow", "darkgreen", "green", "deeppink", "lightblue", "#FFA500", "#00FFFF", "blue", "#698B22")
  
  ggplot(as.data.frame(val_only_transition_types),
         aes(y = Count, axis1 = From, axis2 = To)) +
    geom_alluvium(alpha = 1, aes(fill = From), width = 1/12) +
    
    geom_stratum(stat = "stratum", width = 1/10, aes(fill = after_stat(stratum))) +
    #geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    geom_text(stat = "stratum", aes(label= val_only_transition_types$from_count, position = "dodge")) +
    geom_text(stat = "stratum", aes(label= val_only_transition_types$to_count)) +
    scale_x_discrete(limits = c("From", "To"), expand = c(.05, .05)) +
    
    scale_color_manual(values=cols, aesthetics = c("colour", "fill")) +
    ggtitle("Land cover transition only - validation data") + theme(plot.title = element_text(size=22)) + theme(plot.title=element_text(hjust=0.5))
  theme_void()
  ggsave("Land cover transition only - validation data.pdf", path = "C:\\Master_Thesis\\Output", scale = 2)

  cols <- c("red", "darkgreen", "green", "#00FFFF", "blue", "#698B22", "yellow", "lightblue", "#FFA500")
  
  ggplot(as.data.frame(val_interested_transition_types),
         aes(y = Count, axis1 = From, axis2 = To)) +
    geom_alluvium(alpha = 1, aes(fill = From), width =  1/10, curve_type = "arctangent") +
    
    geom_stratum(stat = "stratum", width = 1/10, aes(fill = after_stat(stratum))) +
    #geom_label(stat = "stratum", aes(label = after_stat(stratum)), how.legend = FALSE, vjust = -0.5) +
    geom_text(stat = "stratum", aes(label= val_interested_transition_types$from_count, position = "dodge")) +
    geom_text(stat = "stratum", aes(label= val_interested_transition_types$to_count)) + #https://stackoverflow.com/questions/60725388/how-to-add-percentage-values-to-strata-in-alluvial-plot-with-ggalluvial
    scale_x_discrete(limits = c("From", "To"), expand = c(.05, .05)) +
    scale_color_manual(values=cols, aesthetics = c("colour", "fill")) +
    ggtitle("Interested land cover transition only - validation data") + theme(plot.title = element_text(size=22)) + theme(plot.title=element_text(hjust=0.5))
  theme_void()
  ggsave("Interested land cover transition only - validation data.pdf", path = "C:\\Master_Thesis\\Output", scale = 2)


