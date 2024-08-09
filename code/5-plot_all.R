# Plot all pops and chromosomes
# part of lidrunner v1.3
# Zane Libke

#library(tidyverse)
library(dplyr)
library(ggplot2)

# Set the directory containing the files
results_dir <- "../results"

# List all files in the directory starting with "LDhat_out_" and ending with ".RData"
filenames <- list.files(path = results_dir, pattern = "^LDhat_out_.*\\.RData$", full.names = TRUE)

print("Plotting data from the following files: ")
print(filenames)

# Load data for each POP_NAME and extract info
all_data <- list()
for (file in filenames) {
  POP_NAME <- sub("^../results/LDhat_out_", "", file)
  POP_NAME <- sub("\\.RData$", "", POP_NAME)
  load(file)
  if (exists("results_all")) {
    results_all$POP_NAME <- POP_NAME
    all_data[[POP_NAME]] <- results_all
    rm(results_all, POP_NAME, file)
    gc()
  } else {
    warning(paste("Variable 'results_all' not found in file:", file))
  }
}

# Combine all data into a single dataframe
combined_df <- bind_rows(all_data)

# Define the desired order of chromosomes
desired_order <- c("2R", "2L", "3R", "3L", "X")

# Convert 'chr' to a factor with specified levels
combined_df$chr <- factor(combined_df$chr, levels = desired_order)

# Check the size of the combined dataframe
print(paste("Combined dataframe size: ", nrow(combined_df), "rows"))

# Parameter string for file names
# this should only pass if all params are equal!
paramstring <- paste0("ww", parameters$window_width, "_spw", parameters$sites_per_window, "_c", parameters$chunk,
                      "_wpc", parameters$windows_per_chunk, "_n", parameters$n, "_r", as.integer(parameters$randsamp))


# Function to calculate average rho per 1MB region for a given dataframe
average_rho_per_region <- function(data, bin_size) {
  data$bin <- floor(data$window_start / bin_size)
  avg_rho <- data %>%
    group_by(chr, bin, POP_NAME) %>%
    summarize(rho = mean(rho, na.rm = TRUE), .groups = 'drop')
  return(avg_rho)
}

# PLOT AVERAGES OVER 1MB
bin_size <- 1000000 
# Calculate average rho per 1MB region for combined dataframe
avg_rho_combined_df <- average_rho_per_region(combined_df, bin_size)

#AVERAGES PLOT
p <- ggplot(avg_rho_combined_df, aes(x = bin, y = rho)) +
  geom_point() +
  facet_grid(POP_NAME ~ chr)+
  labs(x = "Window start (in MB)", y = "avg_rho")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))+
  theme(
    strip.text.x = element_text(size = 6),  # Adjust the facet strip text size for x
    strip.text.y = element_text(size = 6),  # Adjust the facet strip text size for y
    axis.text.x = element_text(size = 6),   # Adjust the x-axis text size
    axis.text.y = element_text(size = 6),    # Adjust the y-axis text size
    panel.spacing = unit(0.1, "lines"),     # Adjust the spacing between facets
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add border to each panel
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"),  # Major gridlines
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey90")  # Minor gridlines
  )

# Increase overall plot size when saving
ggsave("../results/figures/1-avg.png", limitsize = FALSE, plot = p, width = 10, height = 40)  # Adjust width and height as needed

#all windows plot
p <- ggplot(combined_df, aes(x = window_start, y = rho)) +
  geom_point() +
  facet_grid(POP_NAME ~ chr)+
  labs(x = "Window Start (in millions)", y = "Rho")+
  theme_classic()+
  theme(
    strip.text.x = element_text(size = 6),  # Adjust the facet strip text size for x
    strip.text.y = element_text(size = 6),  # Adjust the facet strip text size for y
    axis.text.x = element_text(size = 6),   # Adjust the x-axis text size
    axis.text.y = element_text(size = 6)    # Adjust the y-axis text size
  )

# Increase overall plot size when saving
ggsave("../results/figures/1-all.png", limitsize = FALSE, plot = p, width = 40, height = 80)  # Adjust width and height as needed

