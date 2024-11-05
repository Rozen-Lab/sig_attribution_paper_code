# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())
source("analysis/code/yapsa_analysis.R")
# in_per_sample_cutoff of 0.06 based on Vignette, 
# https://bioconductor.org/packages/release/bioc/vignettes/YAPSA/inst/doc/YAPSA.html
# 3.3.4
run_yapsa("SBS", 0.06)
