# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())
source("analysis/code/QP_analysis.R")
for (thresh in c(0, 0.01, 0.03, 0.06, 0.10)) {
  QP_args = list(trim_less_than_is_fraction = TRUE, trim_less_than = thresh)
  run_QP("SBS", QP_args)
}
