# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())
source("analysis/code/deconstruct_analysis.R")
# 0.06 is the default for deconstructSigs signature.cutoff
run_deconstruct("SBS", 0.06)
