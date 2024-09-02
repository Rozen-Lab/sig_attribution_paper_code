# Run this script with the top level directory as the working directory
source("analysis/code/analysis_utils.R")

time_used <- system.time({
  run_sigpro_syn(
    dataset_name = "SBS",
    python_bin = "/home/gmssgr/miniconda3/envs/sigpro/bin/python",
    run_sigpro_file = "analysis/code/run_sigpro.py",
    seed_in_use = 145879,
    context_type = "96",
    input_root = "analysis/raw_output/SBS/sigpro/syn/input",
    output_root = "analysis/raw_output/SBS/sigpro/syn/output"
  )
})

saveRDS(time_used,
  file = "analysis/raw_output/SBS/sigpro/syn/output/time_used.Rds"
)
