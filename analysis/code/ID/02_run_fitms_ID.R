# Run this script with the top level directory as the working directory
source("analysis/code/fitms_analysis.R")

for (rare_sig_threshold in global_fitms_rare_sig_thresh) {
  run_fitms("ID", rare_sig_threshold)
}
