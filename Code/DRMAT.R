install.packages(c('quantmod','ff','foreign','R.matlab'),dependency=T)
library(tidyverse)
suppressPackageStartupMessages(library(tidyverse))
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(here, devtools, sf, data.table, tidyverse, pbapply, parallel, data.table)
setwd("/home/duyen/Master_Thesis")
source("./Code/Utils/gradual_validation.R")
source("./GitHub_code/BFASTutils.R")
source("./Code/Utils/DRMAT_utils.R")
# add progress bar option to show it in the terminal 
pbo <- pboptions(type="timer")
mycores <- detectCores()

cal_comb_ref <-read.csv("./Intermediate product/cal_compl_data.csv")
cal_ndvi_nonNan <- read.csv("./Data/Data_DRMAT/ndvi_cal.csv")
cal_comb_ref <- tibble::rowid_to_column(cal_comb_ref, "ID")
rm_id <- read.csv("./Data/Data_DRMAT/deleted_sampleid_set.csv")
colnames(rm_id) <- NULL
rm_id <- rm_id[ , 2]
rm_id <- as.list(rm_id)
cal_comb_ref_rm <- cal_comb_ref[!cal_comb_ref$sample_id %in% rm_id,]
bp_BIC_0001 <- DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_BIC_0001.txt", "./Data/Data_DRMAT/bpcm_S_BIC_0001.txt")
bp_BIC_00005 <- DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_BIC_00005.txt", "./Data/Data_DRMAT/bpcm_S_BIC_00005.txt")
bp_BIC_0005 <- DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_BIC_0005.txt", "./Data/Data_DRMAT/bpcm_S_BIC_0005.txt")
bp_BIC_001 <- DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_BIC_001.txt", "./Data/Data_DRMAT/bpcm_S_BIC_001.txt")


bp_AIC_00005 <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_AIC_00005.txt", "./Data/Data_DRMAT/bpcm_S_AIC_00005.txt")
bp_AIC_0005 <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_AIC_0005.txt", "./Data/Data_DRMAT/bpcm_S_AIC_0005.txt")
bp_AIC_0001 <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_AIC_0001.txt", "./Data/Data_DRMAT/bpcm_S_AIC_0001.txt")
bp_AIC_001 <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_AIC_001.txt", "./Data/Data_DRMAT/bpcm_S_AIC_001.txt")
bp_AIC_01 <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_AIC_01.txt", "./Data/Data_DRMAT/bpcm_S_AIC_01.txt")

bp_HQC_00005 <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_HQC_00005.txt", "./Data/Data_DRMAT/bpcm_S_HQC_00005.txt")
bp_HQC_0005 <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_HQC_0005.txt", "./Data/Data_DRMAT/bpcm_S_HQC_0005.txt")
bp_HQC_0001 <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_HQC_0001.txt", "./Data/Data_DRMAT/bpcm_S_HQC_0001.txt")
bp_HQC_001 <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_HQC_001.txt", "./Data/Data_DRMAT/bpcm_S_HQC_001.txt")
bp_HQC_01 <-  DRMAT_conmax("./Data/Data_DRMAT/bpcm_T_HQC_01.txt", "./Data/Data_DRMAT/bpcm_S_HQC_01.txt")
list_of_cal <- rbind( bp_BIC_00005 = bp_BIC_00005,bp_BIC_0001 = bp_BIC_0001, bp_BIC_0005 =  bp_BIC_0005,bp_BIC_001 = bp_BIC_001,
                      bp_AIC_00005 = bp_AIC_00005, bp_AIC_0005 = bp_AIC_0005, bp_AIC_0001 =bp_AIC_0001,
                      bp_AIC_001 = bp_AIC_001, bp_AIC_01 = bp_AIC_01, bp_HQC_00005 = bp_HQC_00005,
                      bp_HQC_0005 = bp_HQC_0005, bp_HQC_0001 =bp_HQC_0001, bp_HQC_001 =bp_HQC_001, bp_HQC_01 =bp_HQC_01)
