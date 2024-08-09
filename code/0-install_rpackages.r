# Install script for lidrunner required R packages
# Zane Libke
# CMEE MSC 2023-4

options(repos = c(CRAN = "https://cran.ma.imperial.ac.uk/"))

packages <- c("ggplot2", "dplyr", "parallel", "hdf5r", "compiler")

install.packages(packages)