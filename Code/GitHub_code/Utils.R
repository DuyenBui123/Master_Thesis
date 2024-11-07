# ---- Scripting preparations ----------------------------------------------------------- 
# load packages 
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(tidyverse, pbapply)

##----------- Time-Series Input File Handling ------------------------------- 

#' function to obtain a data.frame with  
#' NDVI time series and meta data. 
#' observation date columns will be named after those in dfb4. 

#' parameters: 

#' @param dfb4:
#'   data.frame of L8 B4 spectral reflectance data. 
#'   This dataframe should always include the red band of the satellite . 
#'   Columns containing the spectral values should start with
#'   an 'X' and include the date in YYYY-MM-DD format. 


#' @param dfb5:
#'   data.frame of L8 B5 spectral reflectance data. 
#'   This dataframe should always include the NIR band of the satellite.
#'   Columns containing the spectral values should start with
#'   an 'X' and include the date in YYYY-MM-DD format. 

getNDVI <- function(dfb4,dfb5){
  # slice out row 'i', then select all columns starting with an 'X'
  # ensure everything is numeric, and cast the data.frame row to vector
  red <- dfb4 %>% arrange(sample_id) %>%  
    select(starts_with('X')) %>%  
    mutate(across(everything(), as.numeric))
  nir <- dfb5 %>% arrange(sample_id) %>% 
    select(starts_with('X')) %>%  
    mutate(across(everything(), as.numeric))
  
  # get new column names
  newcolnames <- sapply(dfb4 %>% 
                          select(starts_with('X')) %>% 
                          colnames(), 
                        prepareDates)
  
  # calculate the NDVI
  ndvi <- dfb4 %>% 
    # get data.frame with ndvi values
    transmute(round((nir - red)/(nir + red), digits = 4)) %>% 
    # bind all non-sprectral value columns. 
    bind_cols(dfb4 %>% 
                arrange(sample_id) %>%  # make sure to sort by sample_id
                select(!starts_with('X'))) %>% 
    # rearrange the column order
    select(!starts_with('X'), everything()) %>% 
    # place non-numeric columns after the last column 
    relocate(!where(is.numeric), .after = last_col()) %>% 
    # rename colum names that start with an 'X'
    rename_with(.cols = starts_with('X'), ~newcolnames)
  
  return(ndvi)
}


#' Function to standardize Landsat 8 input time series prior to pre-processing. 
#' Individual time-series from the geopackage must be read band by band, 
#' and processed to a standardize ndvi dataframe 
#' with identical columns for every input database (WUR or IIASA).
#' This function takes care of standard meta-data columns, 
#' regrardless of their original database 
#' 
#' @param df: data.frame of raw Landsat 8 band with at least columns
#'            "location_id", "sample_x" and "sample_y"
#' @param refdf: data.frame with reference data. Containing at least columns 
#'            "location_id" and "samle_id"
#' @return data.frame starting with columns "sample_id", 
#'         "centroid_x" and "centroid_y"
standardizeL8TS <- function(df,refdf){
  # make sure id-cols are numeric
  df <- df %>% 
    mutate(across(.cols = ends_with('_id'), 
           ~as.numeric(as.character(.x))))
  refdf <- refdf %>% 
    mutate(across(.cols = ends_with('_id'), 
          ~as.numeric(as.character(.x))))
  
  # prepare join id-columns 
  idjoindf <- refdf %>% 
    distinct(sample_id, location_id)
  
  out <- df %>% 
    # rename columns sample_x and sample_y to centroid_x and centroid_y
    rename_with(~c('centroid_x', 'centroid_y'), 
                .cols = c('sample_x', 'sample_y')) %>% 
    # join the id-columns of refdf to df
    left_join(x=., y= idjoindf, by = "location_id") %>% 
    # put the sample_id column to the front and get rid of the location_id column
    select(sample_id, everything(), -location_id)
  return(out)
}

#' function to apply over a vector to prepare an output vector 
#' with only dates in the format "YYYY-MM-DD"
#' 
#' parameters:
#' @param i: a string assumed to have the form of 
#'    "X2012.12.03_SR_B5" where the date and the band name can differ
prepareDates <- function(i){
  # remove the prefix "X" and split on the underscore "_"
  thedate <- strsplit(str_remove(i, "X"), "_")[[1]][1]
  
  # cast to a date object and format it to "YYYY-MM-DD"
  thedate <- format(as.Date(thedate, "%Y.%m.%d"), "%Y-%m-%d")
  return(thedate)
}

#' Function to obtain a vector with image acquisition dates
#' from the columns of a data frame or an input file. 
#' 
#' @param x data.frame or file name (.csv or .gpkg) of a 
#'          single-layered spatial data file with NDVI values. 
#'          Columns with NDVI observations are assumed 
#'          to start with an "X" and contain their respective dates 
#'          in the column names. 
prepDatesvec <- function(x){
  if(is.character(x)){
    # check if the filename exists
    if(file.exists(x)){
      filename <- x
      # check if the filename is a csv file 
      if(endsWith(filename,".csv")|endsWith(filename,".CSV")){
        l8ndvi <- read.csv(filename)
        
        # check if input filename is a geopackage
      } else if(endsWith(filename,".gpkg")|endsWith(filename,".GPKG")){
        # read the file as a dataframe
        # if more than 1 layer name, pick the layer name that contains 'ndvi'
        lyrnames <- st_layers(filename)$name
        ndvilayer <- ifelse(length(lyrnames) > 1, 
                            lyrnames[grepl('ndvi', lyrnames, ignore.case = TRUE)], 
                            lyrnames)
        l8ndvi <- as.data.frame(sf::st_read(x, layer = ndvilayer, quiet = TRUE))
      } else {stop('Input filename should be either a .csv or a .gpkg')}
    } else {
      stop('please upload the Landsat Time Series Geopackage in the data folder')
    }
  } else if(is.data.frame(x)){
    l8ndvi <- x
  }else{
    stop("x must either be a character filename or a data.frame")
  }
  
  # extract vector with all the date labels
  dateslabels <- l8ndvi %>% colnames() %>% .[str_detect(., pattern = "[X]?[0-2]")]
  
  # extract the dates from the date labels in datesvec
  print("prepare datesvec")
  if(dateslabels %>% startsWith("2") %>% all()){
    datesvec <- as.Date(dateslabels) %>% unname()
  } else {
    datesvec <- pbsapply(dateslabels, FUN = prepareDates, cl=NULL) %>% unname()
    if(datesvec %>% is.na() %>% all()){
      stop('Datesvec is NA')
    }
  }
  return(datesvec)
}

#' Function that preprocesses input file to an analysis ready dataframe
#' 
#' @param inputobject is assumed to be either: 
#'                A character string with a file name corresponding to:
#'                 - a single-layer GPKG file with NDVI values
#'                 - a csv file 
#'                Or a dataframe with columns 
#'                in date-format YYYY-MM-DD, and a tile column 
prepL8NDVIdataframe <- function(inputobject){
  if(is.character(inputobject)){
    filename <- inputobject
    
    # check if the geopackage was uploaded to data folder correctly and read it
    if(file.exists(filename)){
      
      # check if input filename is a csv file: 
      if(endsWith(filename,".csv")|endsWith(filename,".CSV")){
        l8ndvi <- read.csv(filename)
        
        # check if input filename is a geopackage
      } else if(endsWith(filename,".gpkg")|endsWith(filename,".GPKG")){
        # read layer names
        lyrnames <- st_layers(filename)$name
        
        # if more than 1 layer name, pick the layer name that contains 'ndvi'
        ndvilayer <- ifelse(length(lyrnames) > 1, 
                            lyrnames[grepl('ndvi', lyrnames, ignore.case = TRUE)], 
                            lyrnames)
        
        l8ndvi <- sf::st_read(dsn = filename, layer = ndvilayer, quiet = TRUE) %>% 
          as.data.frame()
      } else {
        stop('Filename must be a .csv file or .gpkg file')
      }
    }
  } else if(is.data.frame(inputobject)){
    l8ndvi <- inputobject
  } 
  else{stop("inputobject must either be a character filename or a data.frame")}
  
  # drop the geom column in the dataframe
  if(any(grepl("geo", colnames(l8ndvi)))){
    l8ndvi <- l8ndvi %>% dplyr::select(-contains('geo'))
  }
  
  # get a vector with dates
  datesvec <- l8ndvi %>% prepDatesvec()
  
  # rename columns with NDVI observations to with a name in the form "YYYY-MM-DD"
  l8ndvi <- l8ndvi %>% rename_with(~datesvec, matches("[X]?[0-2]"))
  
  # cast all columns to numeric except for the 'tile' column
  l8ndvi <- l8ndvi %>% mutate(across(-tile, as.numeric))
  
  return(l8ndvi)
}

## ----------- Pre-processing Time Series --------------------------

#' Function to get a freqency table with image acquisitions per year 
#' @param x can either be 
#'          - a character string of the input file (.csv or .gpkg)
#'          - a data.frame with correct date columns
#'          - a vector with dates
getFreqTabByYear <- function(x){
  # prepare a vector with dates
  if(is.character(x)|is.data.frame(x)){
    datesvec <- prepDatesvec(x)
  } else if (is.vector(x)){
    # assume that datesvec was given correctly 
    datesvec <- x
  } else {
    stop("arguments to getFreqTabByYear should either be an input filename or a vector ")
  }
  
  # make look up table with number of observations per year
  tab <- table(cut(as.Date(datesvec, format="%Y-%m-%d"), 'year'))
  # cast it to a well readable data frame 
  tab <- data.frame(Years=format(as.Date(names(tab)), '%Y'),
                    Frequency=as.vector(tab))
  return(tab)
}

## -------------- Output processing -------------------------- 

#' function to join filtered algorithm output to sample id information
#' 
#' function iterates over a vector and takes 
#' columns "sample_id", "centroid_x", "centroid_y" from 
#' a data.frame with the original sample information, 
#' and joins it to the data.frames in the algo_list
#' 
#' parameters:
#' 
#' @param i:          Index of vector to be iterated over
#' @param algo_list:  List with dataframes of selected DBEST output columns 
#' @param dfids:      Dataframe with sample id information.
#'                    Expected columns: "sample_id", "centroid_x", "centroid_y"
joinOutput <- function(i, algo_list, dfids){
  out <- cbind(dfids$sample_id[i], 
               dfids$centroid_x[i], 
               dfids$centroid_y[i], 
               algo_list[[i]])
  if (is.null(colnames(algo_list[[i]])) == TRUE){
    
    colnames(out)[1:4] <- c("sample_id", "centroid_x", "centroid_y", "Rsq")} else {
  colnames(out)[1:3] <- c("sample_id", "centroid_x", "centroid_y")}
  
  return(out)
}

#' Function to calculate the R squared from two vectors
#' formula was obtained from the strucchangeRcpp-package 
#' https://github.com/cran/strucchangeRcpp/blob/2d299537851f77f505cbfc62cfb8418bab07bd20/R/breakpoints.R#L243
#' 
#' @param actualvec: the vector with the actual y values to predict
#' @param predictedvec: the vector with predicted values by the model. 
getR2 <- function(actualvec, predictedvec){
  the_rss <- sum((actualvec - predictedvec) ^ 2)
  mss <- sum((predictedvec - mean(predictedvec))^2)
  the_rsq <- mss/(mss + the_rss)
  return(the_rsq)
}

#' function to attach the name of the algorithm to the output variables in a
#' data.frame.  
#' 
#' @param df: data.frame. If left empty, the function will try to read a df 
#'            from the flename parameter. 
#' @param flename: character. A filename refering to a dataframe. 
#'                 The filename is assumed to have the algorithm name stored
#'                 as the second to last item when the character object 
#'                 is split on underscores and/or dots. Example: "_BFAST.csv"
#'                 When df is not supplied, make sure the flename parameter 
#'                 has the full path to the file. 
#' @param algo_name: character. A characterstring with the algorithm name. 
algonameToVariables <- function(df = "", flename = "", algo_name = ""){
  # if we don't have a data frame, we have to read it first
  if(all(df == "" & flename != "")){
    #' if there was a file name given, 
    #' we read the data frame using the file name
    if(flename != ""){
      df <- read.csv(flename)
    } else{
      # otherwise, raise an error
      stop("No dataframe given, but can't read it from filename either")
    }
  }
  #' if we don't have an algorithm name, 
  #' we have to extract it from the filename
  if(algo_name == ""){
    algo_name <- flename %>%
      # split on the fullstop "." or on the underscore "_"
      strsplit(split = "(\\.|\\_)") %>%
      # vectorize
      unlist() %>%
      # grab the second to last item
      .[length(.)-1]
  }
  
  df <- df %>% rename_with(
    ~case_when(.x == 'value'~ paste0('val', '_',algo_name),
               .x == 'Breakpoint'~ paste0('bp_', algo_name),
               .x == 'Magnitude'~ paste0('mag_', algo_name),
               .x == 'Rsq'~ paste0('r2_', algo_name),
               TRUE ~ .x),
    .cols = everything())
  return(df)
} 


#--------------------- Reference data utilities ---------------------------------

#' Function to read csv files (in this case the ones with reference data)
#' that do not suit the standard read.csv() proceedure.
#' 
#' @param filename: character. The file name to be read. 
#'                  NOTE: as always with file reading, 
#'                  'filename' should contain the full path to the file. 
#' @return df: data.frame with correct and sound classes and column names.
readReferenceCSV <- function(filename){
  # try read the filename in the ordinary way
  df <- read.csv(filename)
  
  # if df has only a single column, we read the file in a different way
  if(ncol(df) == 1){
    # read the file with empty quotation marks
    df <- read.csv(file = filename , sep=",", quote = "", row.names = NULL)
    
    # generate vector with column names
    fileColNames <- readLines(con= filename, n=1) %>% 
      # remove quotation marks by replacing them with nothing
      gsub(pattern = "\"",replacement = "") %>%  
      # split the first line of the file by their commas 
      strsplit(split = ",") %>% 
      # unlist the list of column names
      unlist()
    
    # rename columns 
    colnames(df) <- fileColNames
    
    df <- df %>% 
      # remove all quotation marks in the values of the columns that contain characters 
      mutate(across(where(is.character), 
                    # to remove the quotation marks, we replace them with nothing
                    ~str_replace_all(string =.,  pattern = "\"",replacement =''))) %>% 
      # convert all character columns to their most logical class
      # in this way, character columns with numbers will be set to numeric.
      type.convert(as.is = TRUE)
  }
  return(df)
}

#' Function to add a column with the dominant land cover
#' to the input data.frame (reference data only).
#' 
#' @param df: data.frame with columns for land cover fractions
#' @param lcclasses: character vector. The 'Default' label will set the 
#'                   vector to all known land cover labels in both the 
#'                   WUR and IIASA data sets. Matching of the relevant 
#'                   classes is done automatically. 
#' @return data.frame with character column called "dominant_lc". 
getDominantLC <- function(df, lcclasses = 'Default'){
  # check if there are any land cover classes specified
  if(lcclasses == 'Default'){
    # set the possible land cover classes 
    # "flooded" and "burnt" are classes for temporary events and do not count for 
    # structural land cover change. 
    lcclasses <- c("bare", "crops", "fallow_shifting_cultivation",
                   "grassland", "lichen_and_moss", "not_sure", "shrub", 
                   "snow_and_ice", "tree", "urban_built_up", "water", 
                   "wetland_herbaceous", "fl.grass", "fl.lichen", 
                   "grass", "lichen","snow","trees","urban")
  }
  #check if there was already a dominant_lc column. 
  if("dominant_lc" %in% colnames(df)){
    #' check if the dominant_lc classes 
    #' are all in the list of 'allowed' lc-classes
    #' if that is true, we just return the data.frame 
    if(all(unique(df$dominant_lc) %in% lcclasses)){
      return(df)
    }else{
      df <- df %>% select(-dominant_lc)
    }
  }
  # otherwise we (re)make a 'dominant_lc' column 
  # that indicates which land cover is dominant 
  df$dominant_lc <- df %>% 
    # select all columns with land cover fractions 
    select(dplyr::intersect(colnames(.), lcclasses)) %>% 
    # or select(burnt, flooded, bare, crops, fl.grass, fl.lichen, 
    # grass, lichen, shrub, snow, trees, urban, water)
    {names(.)[max.col(.)]}
  return(df)
}


#' Function to get change columns from reference data.frame
#' @param df: data.frame with land cover fractions or dominant_lc
#' @param idcol: character vector with the name of the ID-column
#' @return data.frame with extra columns "change_at_100m" and "previous_lc"
getLCChangeCols <- function(df, idcol){
  # check if the input df has a column called dominant_lc
  # and whether it is the correct dominant_lc column 
  df <- getDominantLC(df)
  
  # add columns "change_at_100m" and "previous_lc"
  df <- df %>% 
    # group by the ID-column stored in idcol
    group_by(across(all_of(idcol))) %>%    
    # Check if the previous dominant land cover is not equal to the current one. 
    mutate(change_at_100m = lag(dominant_lc, 
                                default = first(dominant_lc)) != dominant_lc) %>% 
    # if there was a change detected, add the previous land cover as a label
    mutate(previous_lc = if_else(
      # if there is a change, store the previous lc label
      # an alternative could be to print pre-post labels
      # paste0(lag(dominant_lc, default = first(dominant_lc)),"-", dominant_lc),
      change_at_100m, # if there is a change
      lag(dominant_lc, default = first(dominant_lc)), # grab the lc class of the previous year
      dominant_lc))  %>% # else, use the current dominant_lc. Alternative: NA_character_
    as.data.frame()
  return(df) 
}

##Handling Benchmark files --------------------------------
#' function to split lines of text file
splitPerfLine <- function(myline){
  if(startsWith(myline, '\t')){
    myline %>% str_replace('\t', '') %>% strsplit(': ')
  }}

#' function to get simple data frame
getSimpleStatDF <- function(i, linelist){
  newcolname = paste0('stat_', i)
  outdf <- data.frame(linelist[[i]][[2]]) %>% 
    rename_with(~newcolname,.cols = everything())
  return(outdf)
}

# function to get data frame from file name
getBenchmarkFromFileName <- function(flename){
  file_lines <- readLines(flename)
  
  linelist <- lapply(file_lines, splitPerfLine) %>% 
    # discard the null-items in the list
    discard(is.null) %>% 
    # unnest list to a clean list of lists
    map(unlist)
  
  #' use getSimpleStatsDF to get a list of simple data frames storing the 
  #' values of the performance statistics under column names 'stat_XX'
  out_df <- lapply(1:length(linelist), getSimpleStatDF, linelist) %>% 
    as.data.frame()
  
  #' create the named list `descr_list`: 
  #' for every colname of df, a description from linelist is defined. 
  descr_list <- linelist %>% 
    lapply(function(i){i[[1]]}) %>% 
    setNames(colnames(out_df))
  
  #' append the descriptions of `descr_list` to the data.frame
  descriptions(out_df) <- descr_list
  
  return(out_df)
}

#' call rbind on a list while retaining the descriptions in the output. 
rbindListWithDescr <- function(mylist){
  out_df <- mylist %>% bind_rows()
  if(!length(descriptions(out_df))){
    descriptions(out_df) <- descriptions(mylist[[1]])
  }
  return(out_df)
}

