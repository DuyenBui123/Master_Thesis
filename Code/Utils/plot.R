# ----  General info --------------------------
#' 
#' Name:             Duyen Bui
#' Study programs:   Master of Geo-information Science 
#' Student number:   1044767
#' 
#' Project:          Thesis: Compare the performance of multivariate and univariate time series models in detecting global land cover change
#' Supervisors:      Dainius Masiliunas & Nandin-Erdene Tsendbazar
#'
#' Script to plot the confusion matrix results of the algorithms
#'
#'_____________________________________________________________________
# install and load package
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(tidyr, ggplot2, gridExtra)
# load external source
debugSource("./Code/Utils/plot_utils.R")
# extract all output files
list_output <- list.files("./Output", pattern = "globe", full.names = TRUE) 
# extract only output files for ndvi val which were interpolated by stlplus
## Plot confusion matrices
list_output_ndvi <- list_output[grepl("NDVI", list_output) & !grepl("SR", list_output)]
cold <- read.csv("./Output/COLD_NDVI_globe.csv" )
# plot with fourfold plots
plot_fourfold_confm(list_output_ndvi)
# plot heatmap
plot_heatmap_confm(list_output_ndvi)
plot_heatmap_confm(list_output_SR)
# plot with barchart
bar_plot <- plot_barchart_confm(list_output_ndvi)
grid.arrange(bar_plot[[1]], bar_plot[[2]], bar_plot[[3]], ncol=3)
dev.off()
# extract only output files for SR val which were interpolated by stlplus
list_output_SR <- list_output[grepl("SR", list_output) & !grepl("NDVI", list_output)]
# plot with fourfold plots
plot_fourfold_confm(list_output_SR)
# plot with barchart
bar_plot_SR <- plot_barchart_confm(list_output_SR)
grid.arrange(bar_plot_SR[[1]], bar_plot_SR[[2]], bar_plot_SR[[3]], ncol=3)
dev.off()


# plot statistics
# read list files for ndvi data
file <- list()
for (file_nr in 1:length(list_output_ndvi)) {
  file_read <- read.csv(list_output_ndvi[file_nr])
  print(file_read)
  file <- append(file, list(file_read))}
# assign the result of each algorithm to an independent variable
Bfastlite <- as.data.frame(file[1])
COLD <- as.data.frame(file[2])
DRMAT <- as.data.frame(file[3])
# create a table that contain the metrices of two algorithms
ndvi <- tibble(
  metrics = c("F1Score", "Precision", "Sensitivity"),
  A_BFASTLite = c(Bfastlite$F1Score, Bfastlite$Precision, Bfastlite$Sensitivity),
  A_COLD = c(COLD$F1Score, COLD$Precision, COLD$Sensitivity),
  A_DRMAT = c(DRMAT$F1Score, DRMAT$Precision, DRMAT$Sensitivity)
  
)
# read list files for SR data
file_SR <- list()
for (file_nr in 1:length(list_output_SR)) {
  file_read <- read.csv(list_output_SR[file_nr])
  print(file_read)
  file_SR <- append(file_SR, list(file_read))}
# assign the result of each algorithm to an independent variable
Bfastlite <- as.data.frame(file_SR[1])
COLD <- as.data.frame(file_SR[2])
DRMAT <- as.data.frame(file_SR[3])
# create a table that contain the metrices of two algorithms
SR <- tibble(
  metrics = c("F1Score", "Precision", "Sensitivity"),
  A_BFASTLite = c(Bfastlite$F1Score, Bfastlite$Precision, Bfastlite$Sensitivity),
  A_COLD = c(COLD$F1Score, COLD$Precision, COLD$Sensitivity),
  A_DRMAT = c(DRMAT$F1Score, DRMAT$Precision, DRMAT$Sensitivity)
  
)
# assign legend, x axis, and y axis
ndvi_plot <- pivot_longer(ndvi, cols = starts_with("A"), names_to = "Algorithms", values_to = "Percentage")
SR_plot <- pivot_longer(SR, cols = starts_with("A"), names_to = "Algorithms", values_to = "Percentage")
# design the size of the plots
options(repr.plot.width = 3, repr.plot.height =4) 
# plot
ndvi_plot1 <- ggplot(ndvi_plot, aes(x = metrics, y = Percentage, fill = Algorithms)) +
  geom_bar(position="dodge", stat="identity", width = 0.75)+ 
  theme_bw() +
  ggtitle("NDVI dataset") + theme(plot.title = element_text(size=40)) + coord_cartesian(ylim = c(0.00,1.00)) +
  theme(plot.title = element_text(hjust = 0.5)) +  theme(legend.position="none")+ theme(text = element_text(size = 34))
ndvi_plot1 + theme(legend.position = "none")
SR_plot1<-ggplot(SR_plot, aes(x = metrics, y = Percentage, fill = Algorithms)) +
  geom_bar(position="dodge", stat="identity", width = 0.75)+ 
  theme_bw() +
  ggtitle("All band dataset") + theme(plot.title = element_text(size=40)) +
  theme(plot.title = element_text(hjust = 0.5)) + theme(text = element_text(size = 34)) + coord_cartesian( ylim = c(0.00,1.00)) +
  theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
        legend.key.height = unit(1.5, 'cm'), #change legend key height
        legend.key.width = unit(1.5, 'cm'), #change legend key width
        legend.title = element_text(size=30), #change legend title font size
        legend.text = element_text(size=30)) #change legend text font size 

# display the plots side by side
grid.arrange(ndvi_plot1, SR_plot1, nrow = 1,  widths = c(1, 1.5))
# save the plots
ggsave("combined_plot.png", plot = grid.arrange(ndvi_plot1, SR_plot1, nrow = 1), width = 5, height = 3)

############## Plot 1 score for ndvi and sr paralelly
plot_parallel_F1(list_output_ndvi, list_output_SR)
# TRY TO ADJUST THE RATIO OF FOURFOLD PLOTS, BUT not yet successed
# matrix_table<- matrix(c(717, 1251, 71795, 797), nrow = 2, byrow = TRUE,
#                       dimnames = list(
#                         Prediction = c("change", "no change"),
#                         Reference = c("change", "no change")
#                       ))
# 
# 
# # Plot with normalized proportions
# fourfoldplot(
#   matrix_table,
#   color = c("#B22222", "#2E8B57"),
#   main = "", conf.level = 0, margin = 1, std = "ind.max"
# )
# 
# # Add annotations
# text(-0.4, 0.4, "TP", cex = 1)
# text(0.4, -0.4, "FP", cex = 1)
# text(0.4, 0.4, "TN", cex = 1)
# text(-0.4, -0.4, "FN", cex = 1)
# dev.off()
# 
# # Load ggplot2
# library(ggplot2)
# 
# # Convert matrix to data frame for plotting
# df <- as.data.frame(as.table(matrix_table))
# colnames(df) <- c("Prediction", "Reference", "Count")
# 
# # Plot using ggplot2
# ggplot(df, aes(x = Reference, y = Count, fill = Prediction)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   scale_fill_manual(values = c("#B22222", "#2E8B57")) +
#   labs(title = "Confusion Matrix Proportions", x = "Reference", y = "Count") +
#   geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.5) +
#   theme_minimal()
# # Apply logarithmic scaling to reduce value disparity
# log_scaled_matrix <- log1p(matrix_table)  # log1p(x) = log(1 + x) to avoid log(0)
# 
# # Plot with log-scaled values
# fourfoldplot(
#   log_scaled_matrix,
#   color = c("#B22222", "#2E8B57"),
#   main = "",
#   conf.level = 0,
#   margin = NULL
# )
# 
# # Add labels for TP, FP, TN, FN
# text(-0.4, 0.4, "TP", cex = 1)
# text(0.4, -0.4, "FP", cex = 1)
# text(0.4, 0.4, "TN", cex = 1)
# text(-0.4, -0.4, "FN", cex = 1)
# 
# # Raw confusion matrix
# matrix_table <- matrix(c(717, 1251, 71795, 797), 
#                        nrow = 2, byrow = TRUE,
#                        dimnames = list(
#                          Prediction = c("change", "no change"),
#                          Reference = c("change", "no change")
#                        ))
# 
# # Normalize the matrix (to avoid a skewed pie chart)
# scaled_matrix <- matrix_table / sum(matrix_table)
# 
# # Plot using fourfoldplot with normalized data
# fourfoldplot(scaled_matrix,
#              color = c("#B22222", "#2E8B57"),
#              main = "",
#              conf.level = 0,
#              margin = NULL)  # No margin normalization
# 
# # Add labels for TP, FP, TN, FN
# text(-0.4, 0.4, "TP", cex = 1)
# text(0.4, -0.4, "FP", cex = 1)
# text(0.4, 0.4, "TN", cex = 1)
# text(-0.4, -0.4, "FN", cex = 1)
# 
# 
# 
# Load necessary libraries
# Load necessary libraries
# Load ggplot2 for visualization

