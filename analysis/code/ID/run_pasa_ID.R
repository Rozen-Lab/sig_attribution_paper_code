# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())
source("analysis/code/pasa_analysis.R")

total_cores <- parallel::detectCores() # 256 on the machine we are using
mc_cores_per_sample <- min(5, total_cores)
pasa_args = list()
pasa_args$mc_cores_per_sample  <- mc_cores_per_sample
pasa_args$num_parallel_samples <- floor(total_cores / mc_cores_per_sample)
pasa_args$seed_in_use = seed = 145879
rm(mc_cores_per_sample, total_cores)
run_pasa(mut_type = "ID",  more_args = pasa_args)
rm(pasa_args)
