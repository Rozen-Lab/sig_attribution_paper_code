# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())
source("analysis/code/deconstruct_analysis.R")
# 0.06 is the default for signature.cutoff
run_deconstruct("DBS", 0.06)
