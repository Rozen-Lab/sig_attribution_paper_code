# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())
source("analysis/code/pasa_analysis.R")

total_cores <- parallel::detectCores() # 256 on the machine we are using
mc_cores_per_sample <- min(5, total_cores)
pasa_SBS_args = list()
pasa_SBS_args$mc_cores_per_sample  <- mc_cores_per_sample
pasa_SBS_args$num_parallel_samples <- floor(total_cores / mc_cores_per_sample)
pasa_SBS_args$seed_in_use = seed = 145879
rm(mc_cores_per_sample, total_cores)
# Note 2024 Sep 9 -- looks like same args as for DBS

run_pasa(mut_type = "SBS",  more_args = pasa_SBS_args)
rm(pasa_SBS_args)
