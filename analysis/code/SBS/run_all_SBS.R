# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())

# This runs everything **except** MSA
 
Sys.setenv(mut_type = "SBS")
Sys.getenv("mut_type")
Sys.setenv(sigpro_context_type = "96")
Sys.getenv("sigpro_context_type")

source("analysis/code/run_all_anything.R")
