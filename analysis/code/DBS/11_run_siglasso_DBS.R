# Run this script with the top level directory as the working directory
stop("sigLASSO does not work on DBS")
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())
source("analysis/code/siglasso_analysis.R")
run_siglasso("DBS", list(use_prior = FALSE))
