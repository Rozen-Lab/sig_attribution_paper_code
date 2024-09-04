# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")

source("analysis/code/deconstruct_analysis.R")

# 0.06 is the default for deconstructSigs
for (cutoff in c(0, 0.03, 0.06, 0.1)) {
  run_deconstruct("ID", cutoff)
}
