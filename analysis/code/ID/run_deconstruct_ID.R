# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(lists = ls())
source("analysis/code/deconstruct_analysis.R")
# 0.06 is the default for deconstructSigs
run_deconstruct("ID", 0.06)
