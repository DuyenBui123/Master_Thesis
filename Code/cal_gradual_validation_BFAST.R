# load packages
# install.packages("pacman"); 
library(pacman)
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(here, devtools, sf, data.table, tidyverse, pbapply, parallel)

# source external functions
here::i_am("Code\\Utils\\gradual_validation.R")
debugSource(here("Code","Utils", "gradual_validation.R"))
# # add progress bar option to show it in the terminal 
pbo <- pboptions(type="timer")
mycores <- detectCores()

NDVI_BIC_025_SandT <-read.csv("C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_SandT_025.csv")
cal_comb_ref <-read.csv("C:\\Master_Thesis\\Intermediate product\\cal_compl_data.csv")

NDVI_BIC_025_SandT <- validateAlgorithmTotal(ref_df = cal_comb_ref, 
                                    algo_df = NDVI_BIC_025_SandT, cl= mycores)
NDVI_BIC_025_SandT_constats <- myFPStats(NDVI_BIC_025_SandT, NewAccuracy = TRUE)

NDVI_BIC_025_T <-read.csv("C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_T_025.csv")
NDVI_BIC_025_T <- validateAlgorithmTotal(ref_df = cal_comb_ref, 
                                             algo_df = NDVI_BIC_025_T, cl= mycores)
NDVI_BIC_025_T_constats <- myFPStats(NDVI_BIC_025_T, NewAccuracy = TRUE)

NDVI_BIC_025_S <-read.csv("C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_S_025.csv")
NDVI_BIC_025_S <- validateAlgorithmTotal(ref_df = cal_comb_ref, 
                                         algo_df = NDVI_BIC_025_S, cl= mycores)
NDVI_BIC_025_S_constats <- myFPStats(NDVI_BIC_025_S, NewAccuracy = TRUE)



NDVI_BIC_1_SandT <-read.csv("C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_SandT_1.csv")
cal_comb_ref <-read.csv("C:\\Master_Thesis\\Intermediate product\\cal_compl_data.csv")

NDVI_BIC_1_SandT <- validateAlgorithmTotal(ref_df = cal_comb_ref, 
                                             algo_df = NDVI_BIC_1_SandT, cl= mycores)
NDVI_BIC_1_SandT_constats <- myFPStats(NDVI_BIC_1_SandT, NewAccuracy = TRUE)

NDVI_BIC_1_T <-read.csv("C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_T_1.csv")
NDVI_BIC_1_T <- validateAlgorithmTotal(ref_df = cal_comb_ref, 
                                         algo_df = NDVI_BIC_1_T, cl= mycores)
NDVI_BIC_1_T_constats <- myFPStats(NDVI_BIC_1_T, NewAccuracy = TRUE)

NDVI_BIC_1_S <-read.csv("C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_BIC_S_1.csv")
NDVI_BIC_1_S <- validateAlgorithmTotal(ref_df = cal_comb_ref, 
                                         algo_df = NDVI_BIC_1_S, cl= mycores)
NDVI_BIC_1_S_constats <- myFPStats(NDVI_BIC_1_S, NewAccuracy = TRUE)
# "LWZ", 1, response ~ harmon + trend
NDVI_LWZ_1_SandT <-read.csv("C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_SandT_1.csv")
cal_comb_ref <-read.csv("C:\\Master_Thesis\\Intermediate product\\cal_compl_data.csv")

NDVI_LWZ_1_SandT <- validateAlgorithmTotal(ref_df = cal_comb_ref, 
                                             algo_df = NDVI_LWZ_1_SandT, cl= mycores)
NDVI_LWZ_1_SandT_constats <- myFPStats(NDVI_LWZ_1_SandT, NewAccuracy = TRUE)

NDVI_LWZ_1_T <-read.csv("C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_T_1.csv")
NDVI_LWZ_1_T <- validateAlgorithmTotal(ref_df = cal_comb_ref, 
                                         algo_df = NDVI_LWZ_1_T, cl= mycores)
NDVI_LWZ_1_T_constats <- myFPStats(NDVI_LWZ_1_T, NewAccuracy = TRUE)

NDVI_LWZ_1_S <-read.csv("C:\\Master_Thesis\\Intermediate product\\cloud_free_product\\_output_BFASTlite_LWZ_S_1.csv")
NDVI_LWZ_1_S <- validateAlgorithmTotal(ref_df = cal_comb_ref, 
                                         algo_df = NDVI_LWZ_1_S, cl= mycores)
NDVI_LWZ_1_S_constats <- myFPStats(NDVI_LWZ_1_S, NewAccuracy = TRUE)


data.table(rbind(
  BF_BIC_025_SandT = myFPStats(NDVI_BIC_025_SandT, NewAccuracy = TRUE), 
 BF_BIC_025_T = myFPStats(NDVI_BIC_025_T, NewAccuracy = TRUE),
 BF_BIC_1_SandT = myFPStats(NDVI_BIC_1_SandT, NewAccuracy = TRUE),
 BF_LWZ_1_SandT = myFPStats(NDVI_LWZ_1_SandT, NewAccuracy = TRUE),
 BF_LWZ_1_T = myFPStats(NDVI_LWZ_1_T, NewAccuracy = TRUE),
 BF_LWZ_1_S = myFPStats(NDVI_LWZ_1_S, NewAccuracy = TRUE)
  
))
