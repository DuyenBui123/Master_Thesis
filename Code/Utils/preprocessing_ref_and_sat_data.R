change_nochange_label <- function(datasource, savepath) {
  ## Define year and number of reference year
  
  datasource$label <- c("no change")
  ## Make a for loop to run to all the unique sample id
  for (sp_id in unique(datasource$sample_id)) {select_sp_id <- datasource[datasource$sample_id == sp_id,] 
  year <- 2015
  ref_nr_year <- 1
  while (ref_nr_year < length(select_sp_id$reference_year)) {
    select_sp_id_1 <- select_sp_id[select_sp_id$reference_year == year,]
    select_sp_id_row1 <- select_sp_id_1 %>% select(bare, crops, grass, lichen, shrub, snow, forest, urban, water, wetland)
    
    select_sp_id_2 <- select_sp_id[select_sp_id$reference_year == year+1,]
    select_sp_id_row2 <- select_sp_id_2 %>% select(bare, crops, grass, lichen, shrub, snow, forest, urban, water, wetland)
    # Calculating change percentage
    if (length(select_sp_id_row1) == length(select_sp_id_row2)) {sub <- abs((select_sp_id_row1 - select_sp_id_row2))
    change_perc <- rowSums(sub, dims = 1)
    print(change_perc)
    ref_nr_year <- ref_nr_year+1
    # Adding threshold
    if(change_perc >= 70) {datasource[datasource$sample_id == sp_id,][ref_nr_year,]$label <- "change"
    } else if (change_perc == 0) {datasource[datasource$sample_id == sp_id,][ref_nr_year,]$label <- "no change"
    } else {datasource[datasource$sample_id == sp_id,][ref_nr_year,]$label <- "small change"}
    year <- year+1
    
    } else{
      print("something wrong with sample id nr:", sp_id)
    }}}
  
  # Write variable into a CSV file
  if (file.exists(savepath)) {
    print("The file exists!")
  } else {
    write.csv(datasource,savepath, row.names = FALSE)
  }
  
}

#' Function to count the number of different types of land cover change transitions
#' @param data is a data.frame that contain all land cover class and their fractional value in percentage 
#' @return return a list of three output containing: all transitions including non-transitions, only transitions and only interested transitions
from_to_transition_type <- function(data) {
  transitionclass_dict <- hash()
  
  for (class in transitionclass) {if (!has.key(class, transitionclass_dict)) {#https://stackoverflow.com/questions/39492564/how-to-create-a-dictionary-and-insert-key-with-a-list-of-value-in-r
    transitionclass_dict[class] <- 0}
  }
  # Assign the class that has the higest percent as the representative class of that row
  for (sample_id in data$sample_id) {sp_id_cal <- data[data$sample_id==sample_id,]
  sel_col_sp_id <- sp_id_cal[,firstclass:lastclass]
  for (row in 1:(length(sp_id_cal$reference_year)-1)) {
    prev_year_class <- colnames(sp_id_cal[row,][,firstclass:lastclass])[apply((sp_id_cal[row,][,firstclass:lastclass]), 1, which.max)]
    consec_year_class <- colnames(sp_id_cal[row+1,][,firstclass:lastclass])[apply((sp_id_cal[row+1,][,firstclass:lastclass]), 1, which.max)]
    from_to_class <- paste0(prev_year_class, "_", consec_year_class)
    transitionclass_dict[[from_to_class]] <- transitionclass_dict[[from_to_class]] + 1
}}
  
  # Extract keys and values from the hash object
  keys <- keys(transitionclass_dict)
  values <- values(transitionclass_dict)
  
  # Create a data frame from the keys and values
  transition_df <- data.frame(
    Transition = keys,
    Count = as.numeric(values)
  )
  transition_df
  # Assign group to the transition class
  transition_df$From <-"none"
  transition_df$To <- "none"
  land_cover_type_trans <- c("bare", "crops", "forest", "grass", "lichen", "shrub", "snow", "urban", "water", "wetland")
  n <- seq(1, 110, by = 10)
  for (i in 1:(length(n)-1)) {
    transition_df$From[n[i]:(n[i+1]-1)] <- land_cover_type_trans[i] # For every 10 row, the group is changed
    transition_df$To[n[i]:(n[i+1]-1)] <- c("bare", "crops", "forest", "grass", "lichen", "shrub", "snow", "urban", "water", "wetland")
    
  }
  # Select only collumns that are necessary. This variable contain all transition and non transition
  all_from_to_class_df <- transition_df %>% select(From, To, Count)
  
  # This variable only contains transition for all classes
  only_transtions <- subset(all_from_to_class_df, all_from_to_class_df[, 1]!=all_from_to_class_df[, 2])
  
  # This variable only contains transition for only interested types
  interested_transition_types <- subset(only_transtions, only_transtions$From %in% c("bare", "crop", "forest", "grass", "urban", "water", "wetland") )
  # Add all three output in a list
  result_lists <- list(all = all_from_to_class_df, only = only_transtions, interested = interested_transition_types )
  return(result_lists)
}

#' Function to count the number of change, no change and small change based on unique sample id
#' @param data is a data.frame that contain all land cover class and a label column of change, no change and small change 

change_nochange_based_sampleid <- function(data) {count_change <- 0
for (sp_id in unique(data$sample_id)) {
  select_sp_id <- data[data$sample_id == sp_id,]
  if((any(select_sp_id$label=="change") == TRUE)) {count_change <-count_change + 1
  
  }}
print(paste0("number of change:", count_change))

count_small_change <- 0
for (sp_id in unique(data$sample_id)) {
  select_sp_id <- data[data$sample_id == sp_id,]
  if((any(select_sp_id$label=="small change") == TRUE)) {count_small_change <-count_small_change + 1
  
  }}
print(paste0("number of small change:", count_small_change))

count_nochange <- 0
for (sp_id in unique(data$sample_id)) {
  select_sp_id <- data[data$sample_id == sp_id,]
  if((all(select_sp_id$label!="change") == TRUE)) {count_nochange <-count_nochange + 1
  
  }}
print(paste0("number of no change:", count_nochange))
}



alluvia_diagram <- function(data, title) {
  # Making alluvia plot
  is_alluvia_form(as.data.frame(data), axes = 1:3, silent = TRUE)
  
  ggplot(as.data.frame(data),
         aes(y = data$Count, axis1 = data$From, axis2 = data$To)) +
    geom_alluvium(aes(fill = data$From), width = 1/12) +
    geom_stratum(width = 1/10, aes(fill = after_stat(stratum)), show.legend = TRUE) +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("From", "To"), expand = c(.05, .05)) +
    
    scale_fill_brewer(type = "qual", palette = "Set1") +
    ggtitle(title)
  theme_void()
}
########################################################## PRE-PROCESSING FOR SAT DATA ##############################################
#' Function to merge two satellite dataset together and 
#' align the satellite data with the reference data to 
#' remove sample ids that are not included in reference data
#'
#' @param data1 is a data.frame of preprocessed WUR satellite dataset with multiple layers that contain at least "location_id" column
#' @param data2 is a data.frame of preprocessed IIASA satellite dataset with multiple layers that contain at least "sample_id", "sample_x", "sample_y" columns
#' @param ref_data is a data.frame of combined referencd data which contains at least "sample_id" and "location_id"
#' @return a data.frame that is a combination of two satellite dataset

# Create a function that merge two satelite dataset together and align the satelite data with the reference data to remove sample ids that are not included in reference data
merge_sat_data <- function(data1, data2, ref_data) {
  data1 <- lapply(data1,function(sat, ref) {sat <- subset(sat, sat$location_id %in% ref$location_id)}, ref_data )
  data2 <- lapply(data2,function(sat, ref) {sat <- subset(sat, sat$sample_id %in% ref$sample_id)}, ref_data )
  for (i in seq_along(data1)) { # Since the data set of IIASA does not have the same format and column order as the data in WUR. Thus a few preprocessing  steps are needed
    data1[[i]]$sample_id <- ref_data$sample_id[match(data1[[i]]$location_id, ref_data$location_id)] # IIASA has sample_id while WUR does not have. Thus ass a sample_id to WUR dataset
  }
  data1 <- lapply(data1, function(df) df[, !names(df) %in% "location_id"]) # Only collect the columns that also exist in the IIASA dataset, so remove location_id column in WUR
  
  rearrage_cols <- function(data) {data %>% select("sample_id", everything())} # Rearrange the columns in WUR to match with the order of columns in IIASA data
  data1 <- lapply(data1, rearrage_cols)

  change_name_sat <- function(data) {colnames(data)[2] <- "centroid_x" # Adjust the column names of sample_x to centroid_x to make it match with WUR column names, similar to ..._y
  colnames(data)[3] <- "centroid_y"
  return(data)}
  data1 <- lapply(data1, change_name_sat)
  
  SR_cal <- Map(rbind, data2, data1) # Finally, merge the two dataset together
  return(SR_cal) }





# # Create a function that merge two satelite dataset together and align the satelite data with the reference data to remove sample ids that are not included in reference data
# merge_sat_data <- function(data1, data2, ref_data) {
#   data1 <- lapply(data1,function(sat, ref) {sat <- subset(sat, sat$location_id %in% ref$location_id)}, ref_data )
#   data2 <- lapply(data2,function(sat, ref) {sat <- subset(sat, sat$sample_id %in% ref$sample_id)}, ref_data )
#   for (i in seq_along(data2)) { # Since the data set of IIASA does not have the same format and column order as the data in WUR. Thus a few preprocessing  steps are needed
#     data2[[i]]$location_id <- ref_data$location_id[match(data2[[i]]$sample_id, ref_data$sample_id)] # IIASA does not have location_id column while the WUR dataset does
#   }
#   data2 <- lapply(data2, function(df) df[, !names(df) %in% "sample_id"]) # Only collect the columns that also exist in the WUR dataset
#   
#   rearrage_cols <- function(data) {data %>% select("location_id", everything())} # Rearrange the columns in IIASA to match with the order of columns in WUR data
#   data2 <- lapply(data2, rearrage_cols)
#   change_name_sat <- function(data) {colnames(data)[2] <- "centroid_x" # Adjust the column names to make it match with WUR column names
#   colnames(data)[3] <- "centroid_y"
#   return(data)}
#   data1 <- lapply(data1, change_name_sat)
#   
#   SR_cal <- Map(rbind, data2, data1) # Finally, merge the two dataset together
#   return(SR_cal) }