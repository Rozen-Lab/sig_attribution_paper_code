# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())
source("analysis/code/yapsa_analysis.R")
# default in_per_sample_cutoff is 0
for (in_per_sample_cutoff in c(0, 0.01, 0.03, 0.06, 0.1)) {
  run_yapsa("DBS", in_per_sample_cutoff)
}

