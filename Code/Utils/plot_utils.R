# install and load package
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(cowplot, ggplot2)

# This function used plot confusion maxtrix results for the algorithms using fourfold plot
#' 
#' @param list_file  a directory path of the output data of all algorithms. they are dataframe and 
#' at least includes TruePositiveCount, FalsePositiveCount, FalseNegativeCount,  TrueNegativeCount 

#' @return do not return any specific variable. The function plots directly the plots
plot_fourfold_confm <- function (list_file) { # create an empty list to collect matrices
  confm <- list()
  for (file.nr in 1:length(list_file)) { # loop through the file list
    # read the file
    file <- read.csv(list_file[file.nr])
    # create a matrix table containing confusion metrics
    matrix_table<- matrix(c(file$TruePositiveCount, file$FalsePositiveCount, file$FalseNegativeCount, file$TrueNegativeCount), nrow = 2, byrow = TRUE,
                          dimnames = list(
                            Prediction = c("change", "no change"),
                            Reference = c("change", "no change")
                          ))
    matrix_table_old <- matrix_table
    matrix_table[1,1] <- round((matrix_table_old[1,1]/sum(matrix_table_old))*100, 2)
    matrix_table[1,2] <- round((matrix_table_old[1,2]/sum(matrix_table_old))*100, 2)
    matrix_table[2,1] <- round((matrix_table_old[2,1]/sum(matrix_table_old))*100, 2)
    matrix_table[2,2] <- round((matrix_table_old[2,2]/sum(matrix_table_old))*100, 2)
    # append the matrix table
    confm <- append(confm,list(matrix_table))
  }
  # create a empty list to collect the names of the algorithms
  algorithmnames <- list()
  for (file.nr in 1:length(list_file)) { # using regex to extract the names
    result <- sub(".*([A-Za-z]+)_.*", "\\1", list_file[file.nr])
    algorithmnames <- append(algorithmnames, result)
  }
  # divide space to plot 
  par(mfrow=c(1,length(list_file)))
  par(cex.main = 1.5, cex.lab = 20, cex.axis = 0.5) # Adjust sizes as needed
  # plot
  for(pltnum in 1:length(confm)) {
    fourfoldplot(confm[[pltnum]],
                 color = c("#B22222", "#2E8B57"),
                 main = "", conf.level = 0, margin = 1, std = "all.max") + 
      text(-0.4,0.4, "TP", cex=2) + 
      text(0.4, -0.4, "TN", cex=2) + 
      text(0.4,0.4, "FP", cex=2) + 
      text(-0.4, -0.4, "FN", cex=2) +
      text (0, 1.5, algorithmnames[pltnum], cex = 2)
    
  }
  
}



# This function used plot confusion maxtrix results for the algorithms using barchart
#' 
#' @param list_file  a directory path of the output data of all algorithms. they are dataframe and 
#' at least includes TruePositiveCount, FalsePositiveCount, FalseNegativeCount,  TrueNegativeCount 

#' @return return all plots in a variable

plot_barchart_confm <- function (list_file) {# create an empty list to collect matrices
  confm <- list()
  for (file in 1:length(list_file)) {# read the file
    file <- read.csv(list_file[file])
    # create a matrix table containing confusion metrics
    matrix_table<- matrix(c(file$TruePositiveCount, file$FalsePositiveCount, file$FalseNegativeCount, file$TrueNegativeCount), nrow = 2, byrow = TRUE)
    
    # Assign labels for TP, TN, FP, FN
    labels <- c("TP", "FN", "FP", "TN")
    # Convert matrices to data frames for plotting with labels
    matrix_table <- as.data.frame(as.table(matrix_table))
    matrix_table_old <- matrix_table
    matrix_table$Freq[1] <-  round((matrix_table_old$Freq[1]/sum(matrix_table_old$Freq))*100, 2)
    matrix_table$Freq[2] <-  round((matrix_table_old$Freq[2]/sum(matrix_table_old$Freq))*100, 2)
    matrix_table$Freq[3] <-  round((matrix_table_old$Freq[3]/sum(matrix_table_old$Freq))*100, 2)
    matrix_table$Freq[4] <-  round((matrix_table_old$Freq[4]/sum(matrix_table_old$Freq))*100, 2)
    # Add TP, TN, FP, FN labels
    matrix_table$Label <- labels
    confm <- append(confm, list(matrix_table))
  }
  
  # Find global maximum value for consistent y-axis
  freq_list <- c()
  for (nr in 1:length(confm)) {
    freq_list <- append(freq_list, confm[[nr]]$Freq)
  }
  global_max <- max(freq_list)
  # create a empty list to collect the names of the algorithms
  algorithmnames <- list()
  for (file in 1:length(list_file)) {# using regex to extract the names
    result <- str_extract(list_file[file], "BFASTlite|COLD|DRMAT|BFASTLite")
    #result <- sub(".*output_([A-Za-z]+)_.*", "\\1", list_file[file])
    algorithmnames <- append(algorithmnames, result)
  }
  plot_list <- list()
  
  for (nr in 1:length(confm)) {
    # Convert to data frame
    df <- as.data.frame(confm[[nr]])
    
    # Create bar plot with labels
    plot <- ggplot(df, aes(x = Label, y = Freq, fill = Label)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_text(aes(label = Freq), 
                position = position_dodge(width = 0.9),  # Align with bars
                vjust = -0.5, 
                size = 10) +  # Adjust text size
      labs(title = algorithmnames[nr], x = "Metrics", y = "Percentage") +
      theme_minimal() +
      scale_fill_manual(values = c("orange", "red", "blue", "green")) +
      scale_y_continuous(limits = c(0, global_max), 
                         breaks = seq(0, global_max, by = round(global_max / 10, 0))) +  # Adjust y-axis tick spacing
      theme(
        plot.title = element_text(size = 40, hjust = 0.5),  # Center title
        legend.position = "none",
        text = element_text(size = 50)
      )
    
    plot_list <- append(plot_list, list(plot))
  }
  return(plot_list)
}

plot_parallel_F1<- function (list_fileNDVI, list_fileSR) {# create an empty list to collect matrices
  ndvi <- c()
  SR <- c()
  names <- c()
  
  for (file in 1:length(list_fileNDVI)) {# read the file
    file_ndvi <- read.csv(list_fileNDVI[file])
    file_sr <- read.csv(list_fileSR[file])
    ndvi <- c(ndvi, file_ndvi$F1)
    SR <- c(SR, file_sr$F1)
    name <- str_extract(list_fileNDVI[file], "BFASTlite|COLD|DRMAT|BFASTLite")
    names <-c(names, name)
  }
  
  value_matrix = matrix(, nrow = 2, ncol = 3)
  value_matrix[1,] = ndvi
  value_matrix[2,] = SR
  
  #Note that the "beside" argument has to be kept "TRUE" in order to place the bars side by side
  barplot(value_matrix, names.arg = names, beside = TRUE, col = c("peachpuff", "skyblue"), legend.text = c("NDVI", "All bands"), ylim = c(0,1), space = c(0,0.5))
  title(main = "F1 scores of NDVI and all band datasets for each algorithm", font.main = 4)
}

library(ggplot2)
library(stringr)
library(cowplot)  # For arranging multiple plots

plot_heatmap_confm <- function(list_file) {
  cm <- list()
  plot_list <- list()
  names <- c()
  
  for (file in list_file) {
    # Read the CSV file
    file_data <- read.csv(file)
    
    # Compute sum of confusion matrix
    sum_cm <- file_data$TruePositiveCount + file_data$FalseNegativeCount + 
      file_data$FalsePositiveCount + file_data$TrueNegativeCount
    
    # Create data frame for confusion matrix
    cm_table <- data.frame(
      Actual = rep(c("Positive", "Negative"), each = 2),
      Predicted = rep(c("Positive", "Negative"), 2),
      Perc = c(round((file_data$TruePositiveCount/sum_cm)*100, 2), 
               round((file_data$FalseNegativeCount/sum_cm)*100, 2), 
               round((file_data$FalsePositiveCount/sum_cm)*100, 2), 
               round((file_data$TrueNegativeCount/sum_cm)*100, 2))  
    )
    
    cm <- append(cm, list(cm_table))
  }
  
  # Extract algorithm names using regex
  for (file.nr in 1:length(list_file)) {
    name <- str_extract(list_file[file.nr], "BFASTlite|COLD|DRMAT|BFASTLite")
    names <- c(names, name)
  }
  
  # Generate confusion matrix plots
  for (nr in 1:length(cm)) {
    plot <- ggplot(data = cm[[nr]], aes(x = Actual, y = Predicted, fill = Perc)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "blue") +
      labs(title = names[nr], x = "Actual", y = "Predicted") +
      geom_text(aes(label = Perc), color = "black", size = 12) +  # Increase number size
      theme_minimal() +
      theme(
        plot.title = element_text(size = 40, hjust = 0.5),  # Centered, large title
        legend.position = "right",
        axis.title = element_text(size = 35),
        axis.text = element_text(size = 30)
      )
    
    plot_list <- append(plot_list, list(plot))
  }
  
  # Arrange multiple plots in a single row
  combined_plot <- plot_grid(plotlist = plot_list, ncol = length(plot_list))
  # Add a title to the whole plot
  final_plot <- ggdraw() + 
    draw_label("Confusion Matrix for Different Algorithms on all band dataset", fontface = 'bold', size = 30,  x = 0.5,
               y = 0.7) +
    draw_plot(combined_plot, y = 0.05, height = 0.5)  # Adjust positioning
  return(final_plot)
}

plot_ts_bp_act <- function(bp_id_uni,NDVI_val, val_comb_ref_remain_change, ts_path ){
  for (id in bp_id_uni) {
    samp_id_bp <- NDVI_val[NDVI_val$sample_id == id,]
    ts_sample_id <- val_comb_ref_remain_change[val_comb_ref_remain_change$sample_id == id,]
    act_bp <- unique(ts_sample_id$change_yr)
    ts_sample_id_total <- val_comb_ref_remain[val_comb_ref_remain$sample_id == id,]
    inter_clfree_ts <- st_read(ts_path)
    ts_samp <- inter_clfree_ts[inter_clfree_ts$sample_id == id,]
    ts_samp <- sf::st_drop_geometry(ts_samp)
    #######
    params <- list(
      preprocessing = list(
        interpolate = FALSE,
        seasonality = ""),
      bfastLiteInputs = list(
        InputTS="",
        scrange=NULL, 
        scsig=0.05, 
        breaks="LWZ",
        sctype="OLS-MOSUM", 
        maginterval=0.1, 
        magcomponent= c("trend"),
        magstat="RMSD", 
        magthreshold=-Inf, 
        coefcomponent=c("trend"),
        coefthresholds=c(0, 0), 
        plot=FALSE, 
        quiet=TRUE, 
        order=3,
        formula=response ~ harmon + trend , 
        TargetYears=NULL,
        seasonfreq=3, 
        breaknumthreshold=Inf, 
        altformula=NULL
      )
    )
    #######
    # load an preprocess the input datafile 
    l8ndvi <- prepL8NDVIdataframe(ts_samp)
    
    # prepare time series ------------------------------------------------------
    ######
    # get frequency table of acquisitions per year
    tab <- getFreqTabByYear(ts_path)
    
    # set the seasonality of the input time series
    if(!is.numeric(params$preprocessing$seasonality)){
      # define the time series frequency by the minimum number of aquisitions in the time series
      params$preprocessing$seasonality <- min(tab$Frequency[2:(length(tab$Frequency)-1)])
    }
    #####
    tslist <- pblapply(1, 
                       FUN = getTrimmedts, 
                       dframe = l8ndvi, 
                       lookuptable = tab, 
                       tsfreq = params$preprocessing$seasonality, 
                       interpolate = params$preprocessing$interpolate, 
                       cl = mycores)
    title_name <- paste0("Time series of sample id:", id, sep=" ")
    title_name <- paste0("Time series of sample id:", id, sep=" ")
    # if (NROW(ts_sample_id) > 0 && all(!is.na(unique(samp_id_bp$Breakpoint))) && all(unique(samp_id_bp$Breakpoint) != c(-9999))) {
    # Plot the first time series (tslist)
    bp <- samp_id_bp$Breakpoint[samp_id_bp$Breakpoint >= 2015 & samp_id_bp$Breakpoint<= 2020]
    par(mar=c(5.1, 4.1, 4.1, 7.1))
    par(mgp = c(2, 1, 0)) 
    my.limits = as.numeric(seq(2013, 2021, by = 1))
    plot(tslist[[1]],                           
         ylim = c(min(as.numeric(tslist[[1]]), na.rm = T), max(as.numeric(tslist[[1]]), na.rm = T)),
         xlab= "",
         ylab="",
         
         type = "l",
         col = "black", # Change to black for better visibility
         lwd = 2       # Line thickness
    )
    title(main = title_name, xlab = "Time series (year)",
          ylab = "NDVI values",
          cex.main = 3,   font.main= 2,
          cex.lab = 1.5, font.lab = 1 )
    axis(1, at = my.limits)
    
    
    
    # Add predicted breakpoints (Red Dashed Lines)
    abline(v = bp, col = "red", lty = 2, lwd = 3, xpd = FALSE) 
    
    # Add actual breakpoints (Blue Dashed Lines)
    abline(v = ts_sample_id$change_yr, col = "blue", lty = 2, lwd = 3, xpd = FALSE) 
    
    # # Add the second time series (non_tslist)
    # lines(non_tslist,                             
    #       type = "l",
    #       col = "green",  # Change to green for distinction
    #       lwd = 2)        # Line thickness
    
    # Add a legend
    legend("topright", inset=c(-0.2,0), cex = 0.8,                       # Position of legend
           legend = c("Non-Inter TS", "Pred bp", "Act bp"), # Labels
           col = c("black", "red", "blue"), # Colors matching plot
           lty = c(1, 1, 2, 2),   # Line types: solid (1) for time series, dashed (2) for breakpoints
           lwd = c(2, 2, 3, 3),   # Line widths
           bty = "n", , xpd = TRUE, 
           title = "Legend")             # No box around the legend
  }
}
# }

# if (NROW(ts_sample_id) > 0 && all(!is.na(unique(samp_id_bp$Breakpoint)))) {
#   # Plot the first time series (tslist)
#   my.limits = as.numeric(seq(2013, 2021, by = 1))
#   plot(tslist[[1]],                           
#        ylim = c(as.numeric(min(tslist[[1]])) -0.1, as.numeric(max(tslist[[1]]))+0.1),
#        
#        type = "l",
#        col = "black", # Change to black for better visibility
#        lwd = 2,       # Line thickness
#        xlab = "Time series (year)",
#        ylab = "NDVI values",
#        main = title_name,
#        cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
#   axis(1, at = my.limits)
#   
#   
#   # Add predicted breakpoints (Red Dashed Lines)
#   abline(v = samp_id_bp$Breakpoint, col = "red", lty = 2, lwd = 3) 
#   
#   # Add actual breakpoints (Blue Dashed Lines)
#   abline(v = ts_sample_id$change_yr, col = "blue", lty = 2, lwd = 3) 
#   
#   
#   # Add a legend
#   legend("topright",                        # Position of legend
#          legend = c("Inter-TS",  "Pred Breakpoint", "Act Breakpoint"), # Labels
#          col = c("black", "red", "blue"), # Colors matching plot
#          lty = c(1,  2, 2),   # Line types: solid (1) for time series, dashed (2) for breakpoints
#          lwd = c(2,  3, 3),   # Line widths
#          bty = "n")             # No box around the legend
# bp <- samp_id_bp$Breakpoint[samp_id_bp$Breakpoint >= 2015 & samp_id_bp$Breakpoint<= 2020]
# # Create a combined data frame for vertical lines
# vline_data <- data.frame(
#   x = c(bp, unique(ts_sample_id$change_yr)),
#   type = c(rep("Pred breakpoint", length(bp)), 
#            rep("Act breakpoint", length(unique(ts_sample_id$change_yr))))
# )
# 
# # Create the plot with a legend
# p <- autoplot(tslist) +
#   geom_vline(data = vline_data, aes(xintercept = x, color = type, linetype = type), size = 1) +
#   scale_x_continuous(breaks = seq(2013, 2022, by = 1)) +
#   labs(
#     title = title_name,
#     x = "Time series (year)",
#     y = "NDVI values",
#     color = "Legend",   # Title of the legend
#     linetype = "Legend" # Ensure linetype also shows in the legend
#   ) +
#   scale_color_manual(values = c("Pred breakpoint" = "red", "Act breakpoint" = "blue")) +
#   scale_linetype_manual(values = c("Pred breakpoint" = "dashed", "Act breakpoint" = "solid")) +
#   theme(plot.title = element_text(hjust = 0.5),
#         text = element_text(size = 25),
#         legend.position = "right") # Adjust legend position if needed
# 
# print(p)
# } else if (NROW(ts_sample_id) > 0 && all(is.na(unique(samp_id_bp$Breakpoint)))){
#   title_name <- paste0("Time series of sample id:", id, "Too much cloud and noise", sep="  ")
#   # Create a combined data frame for vertical lines
#   vline_data <- data.frame(
#     x = c(unique(ts_sample_id$change_yr)),
#     type = c(rep("Act breakpoint", length(unique(ts_sample_id$change_yr))))
#   )
#   
#   # Create the plot with a legend
#   p <- autoplot(tslist) +
#     geom_vline(data = vline_data, aes(xintercept = x, color = type, linetype = type), size = 1) +
#     scale_x_continuous(breaks = seq(2013, 2022, by = 1)) +
#     labs(
#       title = title_name,
#       x = "Time series (year)",
#       y = "NDVI values",
#       color = "Legend",   # Title of the legend
#       linetype = "Legend" # Ensure linetype also shows in the legend
#     ) +
#     scale_color_manual(values = c("Act breakpoint" = "blue")) +
#     scale_linetype_manual(values = c("Act breakpoint" = "solid")) +
#     theme(plot.title = element_text(hjust = 0.5),
#           text = element_text(size = 25),
#           legend.position = "right") # Adjust legend position if needed
#   
#   print(p)
# } else if(NROW(ts_sample_id) >0 && all(unique(samp_id_bp$Breakpoint) == -9999)){
#   title_name <- paste0("Time series of sample id:", id, "no detected breaks", sep="  ")
#   # Create a combined data frame for vertical lines
#   vline_data <- data.frame(
#     x = c(unique(ts_sample_id$change_yr)),
#     type = c(rep("Act breakpoint", length(unique(ts_sample_id$change_yr))))
#   )
#   
#   # Create the plot with a legend
#   p <- autoplot(tslist) +
#     geom_vline(data = vline_data, aes(xintercept = x, color = type, linetype = type), size = 1) +
#     scale_x_continuous(breaks = seq(2013, 2022, by = 1)) +
#     labs(
#       title = title_name,
#       x = "Time series (year)",
#       y = "NDVI values",
#       color = "Legend",   # Title of the legend
#       linetype = "Legend" # Ensure linetype also shows in the legend
#     ) +
#     scale_color_manual(values = c("Act breakpoint" = "blue")) +
#     scale_linetype_manual(values = c("Act breakpoint" = "solid")) +
#     theme(plot.title = element_text(hjust = 0.5),
#           text = element_text(size = 25),
#           legend.position = "right") # Adjust legend position if needed
#   
#   print(p)
# } else {
#   # Create a data frame for vertical lines
#   vline_data <- data.frame(
#     x = c(samp_id_bp$Breakpoint),  # Add more categories if needed
#     type = "Pred breakpoint"  # Label for legend
#   )
#   
#   # Create the plot with a legend
#   p <- autoplot(tslist) +
#     geom_vline(data = vline_data, aes(xintercept = x, color = type, linetype = type), size = 1) +
#     scale_x_continuous(breaks = seq(2013, 2022, by = 1)) +
#     
#     labs(
#       title = title_name,
#       x = "Time series (year)",
#       y = "NDVI values",
#       color = "Legend",   # Legend title
#       linetype = "Legend" # Ensure linetype also appears in legend
#     ) +
#     
#     # Define colors and linetypes in the legend
#     scale_color_manual(values = c("Pred breakpoint" = "red")) +
#     scale_linetype_manual(values = c("Pred breakpoint" = "dashed")) +
#     
#     # Adjust theme and legend position
#     theme(
#       plot.title = element_text(hjust = 0.5),
#       text = element_text(size = 25),
#       legend.position = "right"  # Adjust as needed
#     )
# print(p)
#     }
#   }
#   
# }


plot_ts_bp_act_cold <- function(bp_id_uni,NDVI_val, val_comb_ref_remain_change, ts_path, non_int_ts_path ){
  for (id in bp_id_uni) {
    samp_id_bp <- NDVI_val[NDVI_val$sample_id == id,]
    ts_sample_id <- val_comb_ref_remain_change[val_comb_ref_remain_change$sample_id == id,]
    act_bp <- unique(ts_sample_id$change_yr)
    ts_sample_id_total <- val_comb_ref_remain[val_comb_ref_remain$sample_id == id,]
    inter_clfree_ts <- st_read(ts_path)
    non_clfree_ts <- st_read(non_int_ts_path)
    ts_samp <- inter_clfree_ts[inter_clfree_ts$sample_id == id,]
    non_ts_samp <- non_clfree_ts[non_clfree_ts$sample_id == id,]
    ts_samp <- sf::st_drop_geometry(ts_samp)
    non_ts_samp <- sf::st_drop_geometry(non_ts_samp)
    #######
    
    # Convert column names
    colnames(ts_samp) <- gsub("^X(\\d{4})\\.(\\d{2})\\.(\\d{2})$", "\\1-\\2-\\3", colnames(ts_samp))
    ts_samp <- ts_samp[, !colnames(ts_samp) %in% c("sample_id")  ]
    
    colnames(non_ts_samp) <- gsub("^X(\\d{4})\\.(\\d{2})\\.(\\d{2})$", "\\1-\\2-\\3", colnames(non_ts_samp))
    non_ts_samp <- non_ts_samp[,colnames(non_ts_samp) %in% colnames(ts_samp)]
    
    ts_samp_m <- as.matrix(ts_samp)
    tslist <- ts(ts_samp_m[1,1:186], start = c(2013, 4), end = c(2021, 14),
                 frequency = 22)
    non_ts_samp_m <- as.matrix(non_ts_samp)
    non_tslist <- ts(non_ts_samp_m[1,1:186], start = c(2013, 4), end = c(2021, 14),
                     frequency = 22)
    
    title_name <- paste0("Time series of sample id:", id, sep=" ")
    # if (NROW(ts_sample_id) > 0 && all(!is.na(unique(samp_id_bp$Breakpoint))) && all(unique(samp_id_bp$Breakpoint) != c(-9999))) {
    # Plot the first time series (tslist)
    bp <- samp_id_bp$Breakpoint[samp_id_bp$Breakpoint >= 2015 & samp_id_bp$Breakpoint<= 2020]
    par(mar=c(5.1, 4.1, 4.1, 7.1))
    par(mgp = c(2, 1, 0)) 
    my.limits = as.numeric(seq(2013, 2021, by = 1))
    plot(tslist,                           
         ylim = c(min(as.numeric(tslist)), max(as.numeric(tslist))+0.1),
         xlab= "",
         ylab="",
         
         type = "l",
         col = "black", # Change to black for better visibility
         lwd = 2       # Line thickness
    )
    title(main = title_name, xlab = "Time series (year)",
          ylab = "NDVI values",
          cex.main = 3,   font.main= 2,
          cex.lab = 1.5, font.lab = 1 )
    axis(1, at = my.limits)
    
    
    
    # Add predicted breakpoints (Red Dashed Lines)
    abline(v = bp, col = "red", lty = 2, lwd = 3, xpd = FALSE) 
    
    # Add actual breakpoints (Blue Dashed Lines)
    abline(v = ts_sample_id$change_yr, col = "blue", lty = 2, lwd = 3, xpd = FALSE) 
    
    # Add the second time series (non_tslist)
    lines(non_tslist,                             
          type = "l",
          col = "green",  # Change to green for distinction
          lwd = 2)        # Line thickness
    
    # Add a legend
    legend("topright", inset=c(-0.2,0), cex = 0.8,                       # Position of legend
           legend = c("Inter-TS", "Non-Inter TS", "Pred bp", "Act bp"), # Labels
           col = c("black", "green", "red", "blue"), # Colors matching plot
           lty = c(1, 1, 2, 2),   # Line types: solid (1) for time series, dashed (2) for breakpoints
           lwd = c(2, 2, 3, 3),   # Line widths
           bty = "n", , xpd = TRUE, 
           title = "Legend")             # No box around the legend
  }
  # }
}

plot_ts_bp_act_cold_SR <- function(bp_id_uni,COLD_SRbp, val_comb_ref_remain_change, ts_path, ndvi_path1, noninter_ndvi_path2 ){
  # SR_GPKG = ts_path
  # # Read each name layer of the satellite data
  # SRNames = st_layers(SR_GPKG)$name
  # SR = lapply(SRNames, function(name) st_read(SR_GPKG, layer=name))
  # val_ <- read.csv(ndvi_path1)
  # l8ndvi <- lapply(SR[2:7], function(x) prepL8NDVIdataframe(x))
  # 
  # # prepare time series ------------------------------------------------------
  # 
  # # get frequency table of acquisitions per year
  # tab <- getFreqTabByYear(noninter_ndvi_path2)
  # 
  # # set the seasonality of the input time series
  # if(!is.numeric(params$preprocessing$seasonality)){
  #   # define the time series frequency by the minimum number of aquisitions in the time series
  #   params$preprocessing$seasonality <- min(tab$Frequency[2:(length(tab$Frequency)-1)])
  # }
  # 
  # # make list with timeseries (ts) objects
  # print("Prepare time series list")
  # tslist <- lapply(l8ndvi, function(x) pblapply(1:nrow(x),
  #                                               FUN = getTrimmedts,
  #                                               dframe = x,
  #                                               lookuptable = tab,
  #                                               tsfreq = params$preprocessing$seasonality,
  #                                               interpolate = params$preprocessing$interpolate,
  #                                               cl = mycores))
  for (id in bp_id_uni) {
    index <- which(val_$sample_id ==id)
    samp_id_bp <- COLD_SRbp[COLD_SRbp$sample_id == id,]
    bp <- samp_id_bp$Breakpoint[samp_id_bp$Breakpoint >= 2015 & samp_id_bp$Breakpoint<= 2020]
    val_comb_ref <- val_comb_ref_remain_change[val_comb_ref_remain_change$sample_id == id,]
    
    
    # Add extra space to right of plot area; change clipping to figure
    par(mar=c(5.1, 4.1, 4.1, 7.1))
    par(mgp = c(2, 1, 0)) 
    my.limits = as.numeric(seq(2013, 2021, by = 1))
    title_name <- paste0("Time series of sample id:", id, sep=" ")
    plot(tslist[[1]][[index]],
         ylim = c(min(min(as.integer(tslist[[1]][[index]])),
                      min(as.integer(tslist[[6]][[index]])),
                      min(as.integer(tslist[[5]][[index]]))), max(max(as.integer(tslist[[5]][[index]])),
                                                                  max(as.integer(tslist[[4]][[index]])),
                                                                  max(as.integer(tslist[[6]][[index]])))),
         xlab= "",
         ylab="",
         type = "l",
         col = 1)
    title(main = title_name, xlab = "Time series (year)",
          ylab = "Spectral values",
          cex.main = 3,   font.main= 2,
          cex.lab = 1.7, font.lab = 1 )
    axis(1, at = my.limits)
    
    
    abline(v = bp, col = c("red"),
           lty = c(2), lwd = c(3), xpd = FALSE)
    abline(v = val_comb_ref$change_yr, col = c("blue"),
           lty = c(2), lwd = c(3), xpd = FALSE)
    # Add line graphs of other two dataset    
    lines(tslist[[2]][[index]],
          type = "l",
          col = 2)
    lines(tslist[[3]][[index]],
          type = "l",
          col = 3)
    lines(tslist[[4]][[index]],
          type = "l",
          col = 4)
    lines(tslist[[5]][[index]],
          type = "l",
          col = 7)
    lines(tslist[[6]][[index]],
          type = "l",
          col = 6)
    
    
    # Add a legend
    legend("topright", inset=c(-0.15,0), cex = 0.8, title = "Legend",                     # Position of legend
           legend = c("band 1", "band 2", "band 3", "band 4", "band 5", "band 6" ,"Pred bp", "Act bp"), # Labels
           col = c(1, 2, 3, 4, 7, 6, "red", "blue"), # Colors matching plot
           lty = c(1, 1, 1, 1,1,1, 2, 2),   # Line types: solid (1) for time series, dashed (2) for breakpoints
           lwd = c(2, 2, 2,2,2,2, 3, 3),   # Line widths
           bty = "n", xpd = TRUE)             # No box around the legend
    
    # title_name <- paste0("Time series of sample id:", id, sep=" ")
    # if (NROW(ts_sample_id) > 0 && all(!is.na(unique(samp_id_bp$Breakpoint))) && all(unique(samp_id_bp$Breakpoint) != c(-9999))) {
    #   # Plot the first time series (tslist)
    #   bp <- samp_id_bp$Breakpoint[samp_id_bp$Breakpoint >= 2015 & samp_id_bp$Breakpoint<= 2020]
    #   par(mgp = c(2, 1, 0)) 
    #   my.limits = as.numeric(seq(2013, 2021, by = 1))
    #   plot(tslist,                           
    #        ylim = c(min(as.numeric(tslist)), max(as.numeric(tslist))+0.1),
    #        xlab= "",
    #        ylab="",
    #        
    #        type = "l",
    #        col = "black", # Change to black for better visibility
    #        lwd = 2       # Line thickness
    #   )
    #   title(main = title_name, xlab = "Time series (year)",
    #         ylab = "NDVI values",
    #         cex.main = 3,   font.main= 2,
    #         cex.lab = 1.5, font.lab = 1 )
    #   axis(1, at = my.limits)
    #   
    #   
    #   
    #   # Add predicted breakpoints (Red Dashed Lines)
    #   abline(v = bp, col = "red", lty = 2, lwd = 3) 
    #   
    #   # Add actual breakpoints (Blue Dashed Lines)
    #   abline(v = ts_sample_id$change_yr, col = "blue", lty = 2, lwd = 3) 
    #   
    #   # Add the second time series (non_tslist)
    #   lines(non_tslist,                             
    #         type = "l",
    #         col = "green",  # Change to green for distinction
    #         lwd = 2)        # Line thickness
    #   
    #   # Add a legend
    #   legend("topright",                        # Position of legend
    #          legend = c("Inter-TS", "Non-Inter TS", "Pred Breakpoint", "Act Breakpoint"), # Labels
    #          col = c("black", "green", "red", "blue"), # Colors matching plot
    #          lty = c(1, 1, 2, 2),   # Line types: solid (1) for time series, dashed (2) for breakpoints
    #          lwd = c(2, 2, 3, 3),   # Line widths
    #          bty = "n")             # No box around the legend
  }
}



