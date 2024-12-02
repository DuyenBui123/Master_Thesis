# install and load package
if (!require(cowplot)) install.packages("cowplot")
library(ggplot2)
library(cowplot)
# This function used plot confusion maxtrix results for the algorithms using fourfold plot
#' 
#' @param list_file  a directory path of the output data of all algorithms. they are dataframe and 
#' at least includes TruePositiveCount, FalsePositiveCount, FalseNegativeCount,  TrueNegativeCount 

#' @return do not return any specific variable. The function plots directly the plots
plot_fourfold_confm <- function (list_file) { # create an empty list to collect matrices
  confm <- list()
  for (file in 1:length(list_file)) { # loop through the file list
    # read the file
    file <- read.csv(list_file[file])
    # create a matrix table containing confusion metrics
    matrix_table<- matrix(c(file$TruePositiveCount, file$FalsePositiveCount, file$FalseNegativeCount, file$TrueNegativeCount), nrow = 2, byrow = TRUE,
                          dimnames = list(
                            Prediction = c("change", "no change"),
                            Reference = c("change", "no change")
                          ))
    # append the matrix table
    confm <- append(confm,list(matrix_table))
  }
  # create a empty list to collect the names of the algorithms
  algorithmnames <- list()
  for (file in 1:length(list_file)) { # using regex to extract the names
    result <- sub(".*output_([A-Za-z]+)_.*", "\\1", list_file[file])
    algorithmnames <- append(algorithmnames, result)
  }
  # divide space to plot 
  par(mfrow=c(1,length(list_file)))
  par(cex.main = 1.5, cex.lab = 20, cex.axis = 0.5) # Adjust sizes as needed
  # plot
  for(pltnum in 1:length(confm_per)) {
    fourfoldplot(confm_per[[pltnum]],
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
    result <- sub(".*output_([A-Za-z]+)_.*", "\\1", list_file[file])
    algorithmnames <- append(algorithmnames, result)
  }
  # plot
  plot_list <- list()
  for (nr in 1:length(confm)) {
    # Create barplots for each matrix with shared y-axis limits
    plot <- ggplot(as.data.frame(confm[[nr]]), aes(Label, Freq, fill = Label)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(title = algorithmnames[nr], x = "Metrics", y = "Count") +
      theme_minimal() +
      scale_fill_manual(values = c("orange", "red", "blue", "green")) +
      ylim(0, global_max) + theme(plot.title = element_text(size=40)) +
      theme(plot.title = element_text(hjust = 0.5)) +  theme(legend.position="none")+ theme(text = element_text(size = 50))
    plot_list <- append(plot_list, list(plot))
  }
  return(plot_list)
  
}

