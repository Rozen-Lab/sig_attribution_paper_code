# Run this script with the top level directory as the working directory
source("analysis/code/analysis_utils.R")

time_used <- system.time({
  run_sigpro_syn(
    dataset_name = "ID",
    python_bin = "/home/e0012078/software/miniconda3/bin/python",
    run_sigpro_file = "analysis/code/run_sigpro.py",
    seed_in_use = 145879,
    context_type = "ID",
    input_root = "analysis/raw_output/ID/sigpro_v0.1.7/syn/input",
    output_root = "analysis/raw_output/ID/sigpro_v0.1.7/syn/output"
  )
})

saveRDS(time_used,
        file = "analysis/raw_output/ID/sigpro_v0.1.7/syn/output/time_used.Rds"
)
