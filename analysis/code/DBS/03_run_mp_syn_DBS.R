# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")

source("analysis/code/mp_analysis.R")

run_mp("DBS")
