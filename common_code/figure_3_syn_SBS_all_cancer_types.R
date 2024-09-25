# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())

source("analysis/code/plot_functions.R")

plot_all_cancer_types_merged("SBS", "3")
