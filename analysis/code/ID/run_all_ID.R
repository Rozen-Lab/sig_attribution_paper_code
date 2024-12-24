# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())

 
Sys.setenv(mut_type = "ID")
Sys.getenv("mut_type")
Sys.setenv(sigpro_context_type = "ID")
Sys.getenv("sigpro_context_type")

source("analysis/code/run_all_anything.R")
