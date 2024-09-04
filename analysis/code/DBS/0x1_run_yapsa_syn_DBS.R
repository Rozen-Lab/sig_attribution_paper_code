# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")

source("analysis/code/yapsa_analysis.R")

# default cutoff is 0
for (cutoff in c(0, 0.01, 0.03, 0.06, 0.1)) {
  run_yapsa("DBS", cutoff)
}
