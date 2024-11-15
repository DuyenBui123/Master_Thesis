install.packages("abind")
library(abind)
install.packages("R.matlab")
library(R.matlab)
setwd("/home/duyen/Master_Thesis/")
# setting up the main directory
main_dir <- "/home/duyen/Master_Thesis/Data/"

# setting up the sub directory
sub_dir <- "Data_DRMAT"

# check if sub directory exists 
if (file.exists(sub_dir)){
  
  # specifying the working directory
  setwd(file.path(main_dir, sub_dir))
} else {
  
  # create a new sub directory inside
  # the main path
  dir.create(file.path(main_dir, sub_dir))
  
  # specifying the working directory
  setwd(file.path(main_dir, sub_dir))
}

# Read the CALIBRATION set of satellite data 
SR_GPKG = "./Intermediate product/cloud_free_product/__cloudfree_L8TS__SR_cal.gpkg"
# Read each name layer of the satellite data
SRNames = st_layers(SR_GPKG)$name
SR = lapply(SRNames, function(name) st_read(SR_GPKG, layer=name))


# write dates
outputname_ymd <- paste0(main_dir, "Data_DRMAT/ymd.csv")
date_ts <- colnames(SR_df[[1]])
year_ts <- c()
for (index in 1:length(date_ts)) { year_ts <- append(year_ts,substring(date_ts[index], 2, 5))
}
year_ts_df <- as.data.frame(year_ts)


month_ts <- c()
for (index in 1:length(date_ts)) { month_ts <- append(month_ts,substring(date_ts[index], 7, 8))
}
month_ts_df <- as.data.frame(month_ts)


day_ts <- c()
for (index in 1:length(day_ts)) { day_ts <- append(day_ts,substring(date_ts[index], 10, 11))
}
day_ts_df <- as.data.frame(day_ts)

date_ts_df <- cbind(year_ts_df, month_ts_df, day_ts_df)
colnames(date_ts_df) <- NULL
write.csv(date_ts_df, outputname_ymd)
# Write id
SR_df_id <- lapply(SR,function(x) x <- x[, !names(x) %in% c("tile","centroid_x", "centroid_y")] %>% st_drop_geometry())
outputname <- paste0(main_dir, "Data_DRMAT/id.csv")
sample_id <- as.integer(SR_df_id[[1]][, names(SR_df_id[[1]]) %in% c("sample_id")])
sample_id <- as.data.frame(sample_id)
colnames(sample_id) <- NULL
write.csv(sample_id, outputname)
# Write SR

OutFile = paste0(main_dir, "Data_DRMAT/SR_3D.Rdata")
OutFile_rds = paste0(main_dir, "Data_DRMAT/SR_3D.rds")
OutFile_mat = paste0(main_dir, "Data_DRMAT/SR_3D.mat" )
SR_df <- lapply(SR,function(x) x <- x[, !names(x) %in% c("tile", "sample_id","centroid_x", "centroid_y")] %>% st_drop_geometry())

SR_3D <- abind(SR_df[[2]], SR_df[[3]], SR_df[[4]], SR_df[[5]], SR_df[[6]], SR_df[[7]], along=3)
SR_3D <- unname(SR_3D)

save(SR_3D,file = OutFile)

saveRDS(SR_3D, file = OutFile_rds)
SR_3D_rds <- readRDS(OutFile_rds)
writeMat(OutFile_mat, a = SR_3D_rds)
# write fractional year
date_ts <- sub('X','',date_ts)
date_ts <-  gsub("\\.", "-", date_ts, perl = TRUE)
t <- c()
for (index in 1:length(date_ts)) {
  t <- append(t, sprintf("%.11f",decimal_date(as.POSIXlt(date_ts[index]))))
}
options(digits=11)
t <- as.numeric(t)
t <- as.data.frame(t)
colnames(t) <- NULL
outputname_t <- paste0(main_dir, "/Data_DRMAT/t.csv")
write.csv(t, outputname_t)

# write nivd


