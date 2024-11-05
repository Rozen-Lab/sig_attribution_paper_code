# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())
source("analysis/code/fitms_analysis.R")

for (rare_sig_threshold in global_fitms_rare_sig_thresh) {
  run_fitms("SBS", rare_sig_threshold)
}
