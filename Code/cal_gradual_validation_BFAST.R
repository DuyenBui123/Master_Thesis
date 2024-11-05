# load packages
# install.packages("pacman"); 
library(pacman)
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(here, devtools, sf, data.table, tidyverse, pbapply, parallel, data.table)

# source external functions
here::i_am("Code\\Utils\\gradual_validation.R")
source(here("Code","Utils", "gradual_validation.R"))
# # add progress bar option to show it in the terminal 
pbo <- pboptions(type="timer")
mycores <- detectCores()
setwd("C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\")
NDVI_breakpoint_list <- list.files(path = "C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\", pattern = "^_output_BFASTlite") 

cal_comb_ref <-read.csv("C:\\Master_Thesis\\Intermediate product\\cal_compl_data.csv")
list_breakpoint_output <- data.frame()
for (file in NDVI_breakpoint_list) {
  read_file <- read.csv(file)
  validation <- validateAlgorithmTotal(ref_df = cal_comb_ref, 
                                       algo_df = read_file, cl= mycores)
  stats <- myFPStats(validation, NewAccuracy = TRUE)
  list_breakpoint_output <- rbind(list_breakpoint_output,stats)
}
name_file <- sub(".csv", "", NDVI_breakpoint_list)

list_file <- list()
for (file in name_file) {
  lastchar <- nchar(file)
  file <- substring(file, 18, lastchar)
  list_file <- append(list_file, file)
}

file_name_df <- t(data.frame(list_file, row.names = NULL))
colnames(file_name_df) <- "Parameters"
breakpoint_stats_bfastlite <- cbind(file_name_df,list_breakpoint_output)
write.csv(breakpoint_stats_bfastlite, file = "C:\\Master_Thesis\\Intermediate product\\_stats_breakpoints_bfastlite.csv", row.names = FALSE)