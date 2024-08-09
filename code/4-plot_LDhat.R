#script to take in the .rdata files for each chr produced by 3-run_and_extract_LDhat.R,
#store them as a single df, and plot them
################################### PLOTTING ###################################
#install.packages('ggplot2')
#install.packages('dplyr')
library(ggplot2)
library(dplyr)

# Load in arguments
args <- commandArgs(trailingOnly = TRUE)
POP_NAME <- args[1]

# Load data
filename <- paste("../results/LDhat_out_", POP_NAME, ".RData", sep="")
load(filename)
df <- results_all
paramstring <- paste0("ww", parameters$window_width, "_spw", parameters$sites_per_window, "_c", parameters$chunk,
                      "_wpc", parameters$windows_per_chunk, "_n", parameters$n, "_r", as.integer(parameters$randsamp))

# Get unique chromosomes
chromosomes <- unique(df$chr)

pdf(paste("../results/figures/", POP_NAME, "_", paramstring, ".pdf"))

# PLOTS OF EVERY 100bp WINDOW! NO AVERAGE
for (chromosome in chromosomes) {
  # Filter the dataframe for the current chromosome
  df_chr <- df[df$chr == chromosome, ]
  
  # Create the plot
  p <- plot(df_chr$window_start/1000000, df_chr$rho, type = "p", 
       main = paste("Chromosome", chromosome),
       xlab = "window_start (MB)", ylab = "rho",
       pch = 19, col = "blue")
  print(p)
  # Save the plot to a file
  filename <- paste0("../results/figures/", chromosome, "_",POP_NAME, "_allwind_", paramstring, ".png")
  ggsave(filename, plot = p)
}


# PLOT AVERAGES OVER 1MB
bin_size <- 1000000 

# Function to calculate average rho per 1MB region for a given chromosome
average_rho_per_region <- function(data, bin_size) {
  # Create a new column for the bin
  data$bin <- floor(data$window_start / bin_size)
  
  # Calculate the average rho for each bin
  avg_rho <- aggregate(rho ~ bin, data, mean)
  
  return(avg_rho)
}

# Loop through each unique chromosome and calculate average rho per 1MB region
unique_chromosomes <- unique(df$chr)
for (chromosome in unique_chromosomes) {
  # Filter the dataframe for the current chromosome
  df_chr <- df[df$chr == chromosome, ]
  
  # Calculate average rho per 1MB region
  avg_rho <- average_rho_per_region(df_chr, bin_size)
  
  # Create the plot with points only
  p <- plot(avg_rho$bin, avg_rho$rho, type = "p", 
       main = paste("Avg rho/1MB Chr", chromosome, "_" ,POP_NAME),
       xlab = "window_start(MB)", ylab = "avg rho",
       pch = 19, col = "blue")
  print(p)
  # Save the plot to a file
  filename <- paste0("../results/figures/", chromosome, "_",POP_NAME, "_1MBavg_", paramstring, ".png")
  ggsave(filename, plot = p)
}

dev.off()
