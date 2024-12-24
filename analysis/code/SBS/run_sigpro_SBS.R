# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())
source("analysis/code/sigpro_analysis.R")

sigpro_args = list()
sigpro_args$context_type = "96"
sigpro_args$seed_in_use = global_random_seed
sigpro_args$python_bin = "/home/e0012078/software/miniconda3/bin/python"
run_sigpro("SBS", sigpro_args)
