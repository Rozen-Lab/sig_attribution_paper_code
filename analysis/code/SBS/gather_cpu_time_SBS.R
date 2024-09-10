# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())
source("analysis/code/gather_cpu_time.R")

gather_cpu_info("SBS")
