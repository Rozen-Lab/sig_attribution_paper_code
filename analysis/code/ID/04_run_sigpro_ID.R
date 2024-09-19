# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())
source("analysis/code/sigpro_analysis.R")
sigpro_args = list()
sigpro_args$context_type = "ID"
sigpro_args$seed_in_use = 145879
sigpro_args$python_bin = path.expand("~/miniconda3/envs/sigpro/bin/python")
run_sigpro("ID", sigpro_args)
