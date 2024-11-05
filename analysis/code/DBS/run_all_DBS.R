# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")

Sys.setenv(mut_type = "DBS")
Sys.getenv("mut_type")
Sys.setenv(sigpro_context_type = "DINUC")
Sys.getenv("sigpro_context_type")


source("analysis/code/run_all_anything.R")
