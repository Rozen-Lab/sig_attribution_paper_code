# Run this script with the top level directory as the working directory
source("analysis/code/fitms_analysis.R")

for (rare_sig_threshold in c(0.05, 0.01, 0.03, 0.06)) {
  run_fitms("ID", rare_sig_threshold)
}
