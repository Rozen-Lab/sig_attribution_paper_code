# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())
source("analysis/code/deconstruct_analysis.R")
# 0.06 is the default for signature.cutoff
for (signature.cutoff in c(0, 0.03, 0.06, 0.1)) {
  run_deconstruct("DBS", signature.cutoff)
}

