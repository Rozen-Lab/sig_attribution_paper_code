# Run this script with the top level directory as the working directory
source("analysis/code/analysis_utils.R")

time_used <- system.time({
  run_sigpro_syn(
    dataset_name = "DBS",
    python_bin = "/home/e0012078/software/miniconda3/bin/python",
    run_sigpro_file = "analysis/code/run_sigpro.py",
    seed_in_use = 145879,
    context_type = "DINUC",
    input_root = "analysis/raw_output/DBS/sigpro/syn/input",
    output_root = "analysis/raw_output/DBS/sigpro/syn/output"
  )
})

saveRDS(time_used,
        file = "analysis/raw_output/DBS/sigpro/syn/output/time_used.Rds"
)
